"""Facts about the login host itself, gathered read-only and injectably."""

from __future__ import annotations

import os
import platform
import shutil
import socket
from pathlib import Path
from typing import Callable, Mapping

from .facts import HostFactsV1


def _mem_total_kb(meminfo_path: str = "/proc/meminfo") -> int | None:
    try:
        with open(meminfo_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        return None
    return None


def scratch_candidates(
    env: Mapping[str, str], *, user: str
) -> tuple[str, ...]:
    """Existing directories a scratch could live in, best first."""

    ordered = []
    for candidate in (
        env.get("SCRATCH", ""),
        f"/scratch/{user}",
        "/scratch",
        env.get("TMPDIR", ""),
        "/tmp",
    ):
        candidate = candidate.rstrip("/")
        if candidate and candidate not in ordered:
            ordered.append(candidate)
    return tuple(
        candidate for candidate in ordered if Path(candidate).is_dir()
    )


def gather_host_facts(
    *,
    env: Mapping[str, str] | None = None,
    which: Callable[[str], str | None] = shutil.which,
    meminfo_path: str = "/proc/meminfo",
) -> HostFactsV1:
    env = dict(os.environ) if env is None else env
    user = env.get("USER") or env.get("LOGNAME") or "user"
    return HostFactsV1(
        hostname=socket.gethostname(),
        os_name=platform.system(),
        os_release=platform.release(),
        machine=platform.machine(),
        cpu_count=os.cpu_count(),
        mem_kb=_mem_total_kb(meminfo_path),
        scratch_candidates=scratch_candidates(env, user=user),
        has_module_command=bool(
            which("module") or env.get("MODULESHOME") or env.get("LMOD_DIR")
        ),
    )


__all__ = ["gather_host_facts", "scratch_candidates"]
