"""An approved run handed to a scheduler is still the same run, and the
same one human decision: the job executes the executor's own
continuation inside the allocation and its tail wakes the goal.

Recorded on this host: Slurm accounting storage is disabled, so a
finished job is forgotten after MinJobAge. Nothing here depends on the
scheduler remembering anything -- the run's own record is the evidence.
"""

from __future__ import annotations

import json
import subprocess
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.dispatch import (
    DISPATCH_RECEIPT_FILE,
    EXECUTION_RESULT_FILE,
    DispatchReceiptV1,
    build_dispatch_script,
    dispatch_run_to_scheduler,
    read_dispatch_receipt,
    wait_for_dispatched_run,
)
from chemsmart.agent.driver import DRIVER_FILE, GoalDriver
from chemsmart.settings.server import Server

from .test_the_goal_loop_recovers_or_returns import (
    _engine_stream,
    _envelope_file,
    _planning_session,
    _review_payload,
)


def _slurm_server():
    return Server(
        "canned-slurm",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=6,
        MEM_GB=52,
        NUM_HOURS=24,
        QUEUE_NAME="compute",
        EXTRA_COMMANDS=["ulimit -s unlimited\n"],
    )


def _fake_sbatch(seen):
    def fake_run(argv, **kwargs):
        seen["argv"] = argv
        seen["cwd"] = kwargs.get("cwd")
        return subprocess.CompletedProcess(
            argv, 0, stdout="Submitted batch job 191\n", stderr=""
        )

    return fake_run


def test_the_job_script_runs_the_executor_then_wakes_the_goal(tmp_path):
    server = _slurm_server()
    job = SimpleNamespace(label="goal-g1-cycle-1", PROGRAM=None)
    submitter = server.get_submitter(job)
    script = build_dispatch_script(
        submitter=submitter,
        python="/opt/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    lines = script.splitlines()
    assert lines[0] == "#!/bin/bash"
    assert "#SBATCH --job-name=goal-g1-cycle-1" in lines
    assert "#SBATCH --partition=compute" in lines
    assert "#SBATCH --time=24:00:00" in lines
    assert "ulimit -s unlimited" in lines
    assert "cd $SLURM_SUBMIT_DIR" in lines
    run_line = next(line for line in lines if " agent run " in line)
    assert run_line.startswith("/opt/env/bin/python -m chemsmart agent run")
    assert f"--json > {tmp_path / 'ws' / 'run' / EXECUTION_RESULT_FILE}" in (
        run_line
    )
    wake_line = next(line for line in lines if " agent wake " in line)
    assert wake_line.endswith(f"--workspace {tmp_path / 'ws'} --goal g1")
    assert lines.index(run_line) < lines.index(wake_line)


def test_dispatch_submits_the_script_and_writes_the_receipt(
    tmp_path, monkeypatch
):
    seen: dict = {}
    monkeypatch.setattr(
        "chemsmart.settings.server.subprocess.run", _fake_sbatch(seen)
    )
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.current",
        classmethod(lambda cls: _slurm_server()),
    )
    run_directory = tmp_path / "ws" / ".chemsmart-agent" / "goals" / "g1"
    run_directory = run_directory / "runs" / "cycle-1"
    receipt = dispatch_run_to_scheduler(
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=run_directory,
        goal_id="g1",
        cycle=1,
        python="/opt/env/bin/python",
    )
    assert isinstance(receipt, DispatchReceiptV1)
    assert (receipt.scheduler, receipt.job_id) == ("SLURM", "191")
    assert seen["argv"] == ["sbatch", "chemsmart_sub_goal-g1-cycle-1.sh"]
    assert seen["cwd"] == str(run_directory)
    script = run_directory / "chemsmart_sub_goal-g1-cycle-1.sh"
    assert script.is_file() and " agent wake " in script.read_text()
    stored = read_dispatch_receipt(run_directory)
    assert stored == receipt
    assert json.loads((run_directory / DISPATCH_RECEIPT_FILE).read_text())[
        "wake_command"
    ].endswith("--goal g1")


def test_a_server_without_a_scheduler_refuses_dispatch(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.current",
        classmethod(lambda cls: Server("laptop", SCHEDULER="LOCAL")),
    )
    with pytest.raises(ContractError, match="needs a server profile"):
        dispatch_run_to_scheduler(
            approval_file=tmp_path / "bundle.json",
            workspace=tmp_path,
            run_directory=tmp_path / "run",
            goal_id="g1",
            cycle=1,
        )


def test_waiting_ends_on_the_result_file_or_a_terminal_job(tmp_path):
    run_directory = tmp_path / "run"
    run_directory.mkdir()
    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(
            DispatchReceiptV1(
                scheduler="SLURM",
                job_id="191",
                submitted_at="2026-09-01T00:00:00+00:00",
                submit_command="sbatch x.sh",
                submit_script=str(run_directory / "x.sh"),
                run_directory=str(run_directory),
                approval_file=str(tmp_path / "bundle.json"),
                goal_id="g1",
                cycle=1,
                wake_command="python -m chemsmart agent wake",
            ).public_record()
        ),
        encoding="utf-8",
    )
    states = iter(["RUNNING", "RUNNING", "COMPLETED"])
    slept: list[float] = []

    def runner(argv, **_kw):
        state = next(states)
        return subprocess.CompletedProcess(
            argv, 0, stdout=f"JobId=191 JobState={state} RunTime=00:00:02\n"
        )

    why = wait_for_dispatched_run(
        run_directory, poll_seconds=5, runner=runner, sleep=slept.append
    )
    assert why == "job 191 COMPLETED"
    assert slept == [5, 5]

    (run_directory / EXECUTION_RESULT_FILE).write_text("{}")
    assert (
        wait_for_dispatched_run(run_directory, runner=runner)
        == "result recorded"
    )


def _park_a_goal(tmp_path):
    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)

    def dispatch_run(**kwargs):
        return SimpleNamespace(
            scheduler="SLURM",
            job_id="191",
            submitted_at="2026-09-01T00:00:00+00:00",
            submit_script=str(kwargs["run_directory"] / "sub.sh"),
        )

    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-w1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: _planning_session(
            "live-1", review=_review_payload()
        )(workspace, kw),
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        dispatch_run=dispatch_run,
        dispatch="scheduler",
        server="canned-slurm",
    )
    assert driver.run().settlement == "parked"
    return workspace, driver


def test_the_wake_command_resumes_from_the_recorded_driver(tmp_path):
    from chemsmart.cli.agent import agent

    workspace, driver = _park_a_goal(tmp_path)
    recorded = json.loads((driver.goal_dir / DRIVER_FILE).read_text())
    assert recorded["dispatch"] == "scheduler"
    assert recorded["server"] == "canned-slurm"
    assert recorded["granted_by"] == "claude-owner-delegated-reviewer"

    runner = CliRunner()
    early = runner.invoke(
        agent, ["wake", "--workspace", str(workspace), "--goal", "goal-w1"]
    )
    assert early.exit_code != 0
    assert "has not written its result yet" in early.output

    _engine_stream(tmp_path, driver.run_directory, failed=False)
    (driver.run_directory / EXECUTION_RESULT_FILE).write_text(
        json.dumps({"status": "completed", "analysis_status": ""}),
        encoding="utf-8",
    )
    woken = runner.invoke(
        agent, ["wake", "--workspace", str(workspace), "--goal", "goal-w1"]
    )
    assert woken.exit_code == 0, woken.output
    record = json.loads(woken.output)
    assert record["settlement"] == "achieved"
    assert record["cycles"] == 1

    again = runner.invoke(
        agent, ["wake", "--workspace", str(workspace), "--goal", "goal-w1"]
    )
    assert again.exit_code != 0 and "is settled" in again.output


# --- The scheduler is asked for what the human approved -------------------
#
# The review panel the human reads renders the envelope's resources
# (tui/review.py:_bounds_panel), digest-sealed as resource_sha256. Until
# this section existed, the #SBATCH line came from the operator's server
# profile and the two never met: a live CUHK goal (sm3-run3-cuhk, job
# 2134343) displayed and approved 4 cores / 8 GB and was allocated 32
# cores / 160 GB for a water molecule that ran 31 s at 424 MB peak RSS.
#
# The profile is a ceiling, not a request (owner ruling, 2026-09-16). For
# a sealed job the ceiling is the server's maximum cores and its maximum
# memory less a declared headroom, so a clamped request can never claim a
# node's entire RAM.


def _resources(cores, memory_gb):
    from chemsmart.agent.execution import build_execution_resource_spec

    return build_execution_resource_spec(
        execution_target="run",
        cores=cores,
        memory_gb=memory_gb,
        gpu_count=0,
        scratch_policy="server",
        node_timeout_seconds=900,
    )


def _sbatch_lines(script):
    return [line for line in script.splitlines() if line.startswith("#SBATCH")]


def test_the_scheduler_is_asked_for_the_resources_the_human_approved(
    tmp_path,
):
    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = _slurm_server()  # NUM_CORES=6, MEM_GB=52
    request = resolve_scheduler_request(
        resources=_resources(cores=4, memory_gb=8), server=server, sealed=True
    )
    job = SimpleNamespace(label="goal-g1-cycle-1", PROGRAM=None)
    script = build_dispatch_script(
        submitter=server.get_submitter(job, scheduler_request=request),
        python="/opt/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    lines = _sbatch_lines(script)
    assert "#SBATCH --nodes=1 --ntasks-per-node=4 --mem=8G" in lines, lines
    # The profile still owns queue, wall time and account.
    assert "#SBATCH --partition=compute" in lines
    assert "#SBATCH --time=24:00:00" in lines
    assert not request.clamped


def test_a_sealed_request_is_clamped_to_the_server_s_own_ceiling(tmp_path):
    """The ceiling is max cores and max memory less the declared headroom."""

    from chemsmart.settings.scheduler_request import (
        SEALED_MEMORY_HEADROOM_GB,
        resolve_scheduler_request,
    )

    server = _slurm_server()  # NUM_CORES=6, MEM_GB=52
    request = resolve_scheduler_request(
        resources=_resources(cores=64, memory_gb=300),
        server=server,
        sealed=True,
    )
    assert request.requested_cores == 64
    assert request.requested_memory_gb == pytest.approx(300.0)
    assert request.cores == 6
    assert request.memory_gb == pytest.approx(52.0 - SEALED_MEMORY_HEADROOM_GB)
    assert request.clamped

    job = SimpleNamespace(label="goal-g1-cycle-1", PROGRAM=None)
    script = build_dispatch_script(
        submitter=server.get_submitter(job, scheduler_request=request),
        python="/opt/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    assert "#SBATCH --nodes=1 --ntasks-per-node=6 --mem=46G" in _sbatch_lines(
        script
    )
    # A clamp is a displayed observation naming both numbers, never a
    # silent edit.
    rendered = " | ".join(request.observations)
    assert "cores" in rendered and "64" in rendered and "6" in rendered
    assert "memory" in rendered and "300" in rendered and "46" in rendered


def test_an_unsealed_request_is_the_profile_s_own_numbers(tmp_path):
    """`chemsmart sub` must be byte-identical: no envelope, no headroom."""

    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = _slurm_server()
    request = resolve_scheduler_request(
        resources=None, server=server, sealed=False
    )
    assert (request.cores, request.memory_gb) == (6, 52.0)
    assert not request.clamped

    job = SimpleNamespace(label="j", PROGRAM=None)
    with_request = build_dispatch_script(
        submitter=server.get_submitter(job, scheduler_request=request),
        python="/p",
        approval_file=tmp_path / "b.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    without = build_dispatch_script(
        submitter=server.get_submitter(job),
        python="/p",
        approval_file=tmp_path / "b.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    assert with_request == without


def test_a_ceiling_below_the_headroom_is_refused_not_negated():
    """A 4 GB profile must not resolve to a request for -2 GB."""

    from chemsmart.agent._contracts import ContractError
    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = Server(
        "tiny",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=2,
        MEM_GB=4,
        NUM_HOURS=1,
        QUEUE_NAME="compute",
    )
    with pytest.raises(ContractError, match="headroom"):
        resolve_scheduler_request(
            resources=_resources(cores=1, memory_gb=2),
            server=server,
            sealed=True,
        )


def test_a_clamp_reaches_the_goal_s_own_record_not_only_a_sidecar(
    tmp_path, monkeypatch
):
    """The ledger is what a later process reads.

    A clamp recorded only in dispatch.receipt.json is a clamp the goal's
    own record cannot be audited for -- which is exactly how a live goal
    displayed 4 cores / 8 GB, was allocated 32 / 160, and left nothing in
    its ledger saying so.
    """

    seen: dict = {}
    monkeypatch.setattr(
        "chemsmart.settings.server.subprocess.run", _fake_sbatch(seen)
    )
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.current",
        classmethod(lambda cls: _slurm_server()),
    )
    run_directory = tmp_path / "run"
    receipt = dispatch_run_to_scheduler(
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=run_directory,
        goal_id="g1",
        cycle=1,
        python="/opt/env/bin/python",
        resources=_resources(cores=64, memory_gb=300),
        sealed=True,
    )
    request = receipt.scheduler_request
    assert request is not None, "the receipt records no scheduler request"
    assert request["applied"] == {"cores": 6, "memory_gb": 46, "gpu_count": 0}
    assert request["requested"]["cores"] == 64
    assert request["ceiling"]["memory_gb"] == 46
    assert request["clamped"] is True
    assert len(request["observations"]) == 2

    # The driver copies a subset of the receipt into the ledger row; the
    # scheduler request must be inside that subset, so drive the filter
    # rather than trusting it.
    from chemsmart.agent.driver import _record_of

    payload_keys = {
        "scheduler",
        "job_id",
        "submitted_at",
        "submit_script",
        "wake_job_id",
        "scheduler_request",
    }
    carried = {
        key: value
        for key, value in _record_of(receipt).items()
        if key in payload_keys
    }
    assert "scheduler_request" in carried, (
        "the run_dispatched ledger row drops the scheduler request, so a "
        "clamp is invisible to every later reader of the goal record"
    )
    assert carried["scheduler_request"]["clamped"] is True


def test_the_resolver_is_the_only_author_of_an_sbatch_resource_line():
    """No path may still read the profile's cores or memory directly.

    A resolver that one writer bypasses is a resolver that is true in the
    tests and false on whichever path was missed.
    """

    import re
    from pathlib import Path

    source = Path("chemsmart/settings/submitters.py").read_text()
    offenders = [
        line.strip()
        for line in source.splitlines()
        if re.search(r"#SBATCH|#PBS|#BSUB|#PJM", line)
        and re.search(r"server\.(num_cores|mem_gb|num_gpus)", line)
    ]
    assert not offenders, (
        "these scheduler directives still author resources from the "
        f"profile instead of the resolved request: {offenders}"
    )
