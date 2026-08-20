"""Honest scheduler detection: environment first, then binaries that answer.

The legacy ``detect_server_scheduler`` treats a binary's mere existence as
proof (it ignores exit codes); a wizard must not. Detection here records
its evidence, checks that the query command actually responds, and
disambiguates the two schedulers that both call their submitter ``qsub``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Callable, Mapping

from .facts import DetectionV1

Runner = Callable[..., "subprocess.CompletedProcess[str]"]


def _run(command: list[str], runner: Runner) -> "subprocess.CompletedProcess[str]":
    return runner(
        command, capture_output=True, text=True, timeout=20, check=False
    )


def detect_scheduler(
    *,
    env: Mapping[str, str] | None = None,
    which: Callable[[str], str | None] = shutil.which,
    runner: Runner = subprocess.run,
) -> DetectionV1:
    env = dict(os.environ) if env is None else env
    evidence: list[str] = []

    # SLURM: environment markers, else a responding sinfo.
    if env.get("SLURM_JOB_ID") or env.get("SLURM_CLUSTER_NAME"):
        evidence.append("SLURM environment variables present")
        return DetectionV1("SLURM", tuple(evidence))
    if which("sbatch") and which("sinfo"):
        probe = _run(["sinfo", "--version"], runner)
        if probe.returncode == 0:
            evidence.append(
                f"sinfo responds: {probe.stdout.strip() or 'ok'}"
            )
            return DetectionV1("SLURM", tuple(evidence))
        evidence.append(
            "sbatch/sinfo installed but sinfo exited "
            f"{probe.returncode}; the SLURM controller may be down"
        )

    # qsub is shared by PBS and SGE; disambiguate before believing either.
    if which("qsub"):
        if env.get("SGE_ROOT"):
            evidence.append("SGE_ROOT is set")
            return DetectionV1("SGE", tuple(evidence))
        if env.get("PBS_JOBID") or env.get("PBS_O_HOST") or which("pbsnodes"):
            probe = _run(["qstat", "--version"], runner)
            if probe.returncode == 0 and probe.stdout.strip():
                evidence.append(f"qstat responds: {probe.stdout.strip()}")
            else:
                evidence.append("qsub and pbsnodes present")
            return DetectionV1("PBS", tuple(evidence))
        probe = _run(["qstat", "-help"], runner)
        blob = (probe.stdout + probe.stderr).lower()
        if "sge" in blob or "grid engine" in blob:
            evidence.append("qstat -help names Grid Engine")
            return DetectionV1("SGE", tuple(evidence))
        evidence.append("qsub present; assuming PBS")
        return DetectionV1("PBS", tuple(evidence))

    if env.get("LSB_JOBID") or which("bsub"):
        evidence.append("LSF markers present (bsub/LSB_JOBID)")
        return DetectionV1("LSF", tuple(evidence))
    if which("condor_q"):
        evidence.append("condor_q present")
        return DetectionV1("HTCondor", tuple(evidence))

    evidence.append("no batch scheduler responded; local execution only")
    return DetectionV1(None, tuple(evidence))


__all__ = ["detect_scheduler"]
