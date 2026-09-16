"""What the scheduler allocated is what the engine is told it may use.

The server profile decides the allocation, and the ``#SBATCH`` line was
taught to read it while the executor went on writing
``execution-server.yaml`` from ``bundle.execution_resources`` -- the
episode's own numbers. So a job allocated 64 cores told ORCA to use 4, or
a job allocated 46 GB told it 300: one question, two answers, and the
engine believing the wrong one.

That is one question with two answers, which is the defect class the round
exists to remove, introduced by the commit that was removing it. It is the
worst shape of it, too: the failure surfaces as an OOM kill or a native
program error, lands in a repairable terminal state, and offers the
session a *scientific* repair menu for what was a host arithmetic
decision.

Found by an adversarial review of the implementation (checkpoint B,
2026-09-16), not by the suite, which asserted only what the script said.
"""

from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.execution import build_execution_resource_spec
from chemsmart.agent.live_session import _write_execution_server_profile


def _resources(cores: int, memory_gb: float):
    return build_execution_resource_spec(
        execution_target="run",
        cores=cores,
        memory_gb=memory_gb,
        gpu_count=0,
        scratch_policy="server",
        node_timeout_seconds=900,
    )


def _profile_values(text: str) -> dict[str, str]:
    values = {}
    for line in text.splitlines():
        stripped = line.strip()
        for key in ("NUM_CORES", "NUM_THREADS", "MEM_GB"):
            if stripped.startswith(f"{key}:"):
                values[key] = stripped.split(":", 1)[1].strip()
    return values


def _write_receipt(run_directory: Path, *, applied, requested, ceiling):
    from chemsmart.agent.dispatch import DISPATCH_RECEIPT_FILE

    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(
            {
                "scheduler": "SLURM",
                "job_id": "191",
                "submitted_at": "2026-09-16T00:00:00+00:00",
                "submit_command": "sbatch x.sh",
                "submit_script": str(run_directory / "x.sh"),
                "run_directory": str(run_directory),
                "approval_file": str(run_directory / "bundle.json"),
                "goal_id": "g1",
                "cycle": 1,
                "wake_command": "python -m chemsmart agent wake",
                "schema_version": "chemsmart.goal-dispatch-receipt.v1",
                "scheduler_request": {
                    "applied": applied,
                    "requested": requested,
                    "ceiling": ceiling,
                    "sealed": True,
                    "clamped": applied != requested,
                    "observations": ["cores clamped", "memory clamped"],
                },
            }
        ),
        encoding="utf-8",
    )


def test_the_granted_allocation_is_what_the_engine_is_told(tmp_path):
    """One set of numbers for the allocation and the program inside it."""

    run_directory = tmp_path / "run"
    run_directory.mkdir()
    _write_receipt(
        run_directory,
        applied={"cores": 32, "memory_gb": 154, "gpu_count": 0},
        requested={"cores": 64, "memory_gb": 300, "gpu_count": 0},
        ceiling={"cores": 32, "memory_gb": 154, "gpu_count": 0},
    )

    profile = _write_execution_server_profile(
        run_directory, _resources(cores=64, memory_gb=300)
    )
    values = _profile_values(profile.read_text())
    assert values["NUM_CORES"] == "32", (
        "the engine was told a core count the allocation does not have: "
        f"{values}"
    )
    assert values["NUM_THREADS"] == "32"
    assert values["MEM_GB"] == "154", (
        "the engine was told a memory the cgroup does not allow, which is "
        f"an OOM kill wearing a scientific failure's word: {values}"
    )


def test_an_unclamped_dispatch_leaves_the_engine_exactly_as_approved(
    tmp_path,
):
    run_directory = tmp_path / "run"
    run_directory.mkdir()
    _write_receipt(
        run_directory,
        applied={"cores": 4, "memory_gb": 8, "gpu_count": 0},
        requested={"cores": 4, "memory_gb": 8, "gpu_count": 0},
        ceiling={"cores": 32, "memory_gb": 154, "gpu_count": 0},
    )
    profile = _write_execution_server_profile(
        run_directory, _resources(cores=4, memory_gb=8)
    )
    values = _profile_values(profile.read_text())
    assert (values["NUM_CORES"], values["MEM_GB"]) == ("4", "8")


def test_a_local_run_has_no_receipt_and_is_unchanged(tmp_path):
    """No scheduler, no clamp: the envelope is the whole truth."""

    run_directory = tmp_path / "run"
    run_directory.mkdir()
    profile = _write_execution_server_profile(
        run_directory, _resources(cores=4, memory_gb=8)
    )
    values = _profile_values(profile.read_text())
    assert (values["NUM_CORES"], values["NUM_THREADS"], values["MEM_GB"]) == (
        "4",
        "4",
        "8",
    )


def test_a_larger_allocation_than_the_episode_asked_for_is_used(tmp_path):
    """The profile is the authority in both directions.

    A server profile of 64 cores for an episode that asked for 4 means the
    job has 64, so the engine is told 64 -- otherwise 60 allocated cores
    sit idle while the operator's own setting is ignored.
    """

    run_directory = tmp_path / "run"
    run_directory.mkdir()
    _write_receipt(
        run_directory,
        applied={"cores": 64, "memory_gb": 154, "gpu_count": 0},
        requested={"cores": 4, "memory_gb": 8, "gpu_count": 0},
        ceiling={"cores": 64, "memory_gb": 154, "gpu_count": 0},
    )
    profile = _write_execution_server_profile(
        run_directory, _resources(cores=4, memory_gb=8)
    )
    values = _profile_values(profile.read_text())
    assert (values["NUM_CORES"], values["NUM_THREADS"], values["MEM_GB"]) == (
        "64",
        "64",
        "154",
    )
