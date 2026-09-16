"""Two elements of one wave must derive the same allocation.

`_allocated_execution_resources` reads `dispatch.receipt.json` -- the
host's own record of what the scheduler granted -- and that receipt was
written *after* both submissions. The manifest is written first precisely
because elements start promptly, and on CUHK four of them started in the
same second.

So an element inside that window finds no receipt, falls back to the
*approval's* resources, and publishes those; an element after it derives
the granted allocation and finds bytes that disagree. Making the shared
write atomic changed the symptom and not the cause: instead of three
elements dying of `FileExistsError`, the one element that runs is the one
running under the wrong numbers -- `%pal nprocs 64 end` inside a 32-core
cgroup, arriving as a repairable *scientific* terminal state for a host
ordering mistake. That is the defect `8953ef43` closed, re-entered
through the cohort path.

The allocation is known before any `sbatch`, so it is written before any
`sbatch`, for the same reason the manifest is.
"""

from __future__ import annotations

import json
import subprocess
from types import SimpleNamespace

import pytest

from chemsmart.agent.dispatch import (
    DISPATCH_RECEIPT_FILE,
    dispatch_run_to_scheduler,
)
from chemsmart.agent.execution import build_execution_resource_spec
from chemsmart.settings.server import Server

_NODES = ("a1", "a2", "a3")
_BUNDLE = "7" * 64


def _server():
    return Server(
        "canned-slurm",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=16,
        MEM_GB=40,
        NUM_HOURS=4,
        QUEUE_NAME="compute",
    )


@pytest.fixture(autouse=True)
def _profile(monkeypatch):
    monkeypatch.setattr(
        Server, "from_servername", classmethod(lambda cls, name: _server())
    )


def test_an_element_that_starts_before_the_job_id_lands_reads_the_grant(
    tmp_path, monkeypatch
):
    """The window is the submission itself."""

    from chemsmart.agent.live_session import (
        _allocated_execution_resources,
    )

    run_directory = tmp_path / "run"
    bundle = tmp_path / "bundle.json"
    bundle.write_text(json.dumps({"bundle_sha256": _BUNDLE}))
    # The approval asked for far more than this profile grants, which is
    # ordinary: the profile is the authority and nothing is refused.
    approved = build_execution_resource_spec(
        execution_target="run",
        cores=64,
        memory_gb=300,
        gpu_count=0,
        scratch_policy="server",
        node_timeout_seconds=900,
    )

    seen: list[SimpleNamespace] = []

    def fake_run(argv, **kwargs):
        # Exactly what a promptly-started element sees: the job exists,
        # and this process has not finished writing its receipt.
        seen.append(
            _allocated_execution_resources(run_directory, approved)
        )
        return subprocess.CompletedProcess(
            argv, 0, stdout=f"Submitted batch job 4{len(seen)}2\n", stderr=""
        )

    monkeypatch.setattr(subprocess, "run", fake_run)
    dispatch_run_to_scheduler(
        approval_file=bundle,
        workspace=tmp_path / "ws",
        run_directory=run_directory,
        goal_id="g1",
        cycle=1,
        server="canned-slurm",
        python="/opt/env/bin/python",
        resources=approved,
        cohort_node_ids=_NODES,
    )

    assert seen, "the submission never happened"
    granted = _allocated_execution_resources(run_directory, approved)
    for index, during in enumerate(seen):
        assert during.cores == granted.cores, (
            f"an element starting during submission {index} would derive "
            f"{during.cores} cores while the wave was granted "
            f"{granted.cores}: one element runs outside its allocation "
            "and the rest die of conflicting bytes"
        )
        assert during.memory_gb == granted.memory_gb

    # And the grant is the profile's, not the approval's.
    assert granted.cores == 16
    assert granted.cores != approved.cores


def test_the_receipt_still_names_the_jobs_when_it_is_complete(tmp_path):
    """Writing early must not lose what is only known late."""

    import subprocess as _sp

    run_directory = tmp_path / "run"
    bundle = tmp_path / "bundle.json"
    bundle.write_text(json.dumps({"bundle_sha256": _BUNDLE}))
    calls = {"n": 0}

    def fake_run(argv, **kwargs):
        calls["n"] += 1
        return _sp.CompletedProcess(
            argv,
            0,
            stdout=f"Submitted batch job 4{calls['n']}2\n",
            stderr="",
        )

    original = _sp.run
    _sp.run = fake_run
    try:
        receipt = dispatch_run_to_scheduler(
            approval_file=bundle,
            workspace=tmp_path / "ws",
            run_directory=run_directory,
            goal_id="g1",
            cycle=1,
            server="canned-slurm",
            python="/opt/env/bin/python",
            cohort_node_ids=_NODES,
        )
    finally:
        _sp.run = original

    record = json.loads(
        (run_directory / DISPATCH_RECEIPT_FILE).read_text(encoding="utf-8")
    )
    assert record["job_id"] == "412" == receipt.job_id
    assert record["wake_job_id"] == "422"
    assert record["cohort_node_ids"] == list(_NODES)
