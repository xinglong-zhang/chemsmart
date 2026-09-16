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
    # The profile still owns queue and account.
    assert "#SBATCH --partition=compute" in lines
    # Wall time is the envelope's too. With no envelope record the node
    # timeout is the only bound the host has -- 900 s, so one hour -- and
    # it is used rather than the operator's 24 h default. A 31-second
    # water job holding a day of wall clock is not backfill-eligible and
    # charges fairshare against a shared reservation for the whole day.
    assert "#SBATCH --time=1:00:00" in lines
    assert not request.clamped


def test_a_sealed_request_above_the_ceiling_is_refused_not_shrunk(
    tmp_path,
):
    """The host does not choose a method.

    How much memory a correlated calculation gets selects the algorithm --
    integral-direct against conventional, disk against in-core -- and
    sometimes decides whether it runs at all. Shrinking it silently makes
    that choice on the scientist's behalf and delivers the consequence as
    an OOM kill wearing a native-program failure's word, in a repairable
    terminal state, with a scientific repair menu offered for a host
    arithmetic decision (owner ruling, 2026-09-16, revised on evidence).
    """

    from chemsmart.agent._contracts import ContractError
    from chemsmart.settings.scheduler_request import (
        SEALED_MEMORY_HEADROOM_GB,
        resolve_scheduler_request,
    )

    server = _slurm_server()  # NUM_CORES=6, MEM_GB=52
    with pytest.raises(ContractError) as excinfo:
        resolve_scheduler_request(
            resources=_resources(cores=64, memory_gb=300),
            server=server,
            sealed=True,
        )
    message = str(excinfo.value)
    # Both numbers, in both dimensions, and the headroom that made the
    # memory ceiling what it is.
    assert "64" in message and "6" in message
    assert "300" in message
    assert str(int(52 - SEALED_MEMORY_HEADROOM_GB)) in message
    assert "headroom" in message
    # And the routes onward, because a refusal that teaches nothing is a
    # refusal a session meets five times in a row.
    assert "envelope" in message and "server profile" in message


def test_a_request_at_the_ceiling_is_allowed(tmp_path):
    from chemsmart.settings.scheduler_request import (
        SEALED_MEMORY_HEADROOM_GB,
        resolve_scheduler_request,
    )

    server = _slurm_server()
    request = resolve_scheduler_request(
        resources=_resources(
            cores=6, memory_gb=52 - SEALED_MEMORY_HEADROOM_GB
        ),
        server=server,
        sealed=True,
    )
    assert (request.cores, request.memory_gb) == (6, 46)
    assert not request.clamped


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
        resources=_resources(cores=4, memory_gb=8),
        sealed=True,
    )
    request = receipt.scheduler_request
    assert request is not None, "the receipt records no scheduler request"
    assert request["applied"]["cores"] == 4
    assert request["applied"]["memory_gb"] == 8
    assert request["ceiling"]["memory_gb"] == 46
    assert request["ceiling"]["cores"] == 6

    # The driver copies a subset of the receipt into the ledger row; the
    # scheduler request must be inside that subset, so drive the filter
    # rather than trusting it.
    from chemsmart.agent.driver import _dispatch_ledger_keys, _record_of

    # The filter is read from the driver rather than restated here: a
    # test that re-types the producer's key set is the pattern that let
    # the array rot, and one of the driver's keys ("wake_job_id") has no
    # producer at all -- asserting it would pin a dead literal as if it
    # were a contract.
    carried = {
        key: value
        for key, value in _record_of(receipt).items()
        if key in _dispatch_ledger_keys()
    }
    assert "scheduler_request" in carried, (
        "the run_dispatched ledger row drops the scheduler request, so a "
        "clamp is invisible to every later reader of the goal record"
    )
    assert carried["scheduler_request"]["ceiling"]["memory_gb"] == 46
    # wake_command is the only durable record of how this goal is meant
    # to be resumed, and the filter used to drop it.
    assert carried["wake_command"].endswith("--goal g1")


def test_the_resolver_is_the_only_author_of_a_scheduler_resource_line():
    """No submitter may read a resource off the profile directly.

    The first form of this guard matched a directive prefix and a
    ``server.<attr>`` on the *same physical line*, from a relative path.
    It missed three things at once: the PBS block spells its request over
    two lines, so the line carrying ``mem=`` has no ``#PBS`` token;
    ``num_nodes`` and ``num_threads`` were not in the alternation, and
    ``num_nodes`` was a live AttributeError; and reading
    ``Path("chemsmart/...")`` made the whole check depend on the working
    directory. A lint that can be walked around is a lint that will be.

    So this reads the module by AST, from the package itself, and refuses
    any attribute access of the form ``self.server.<resource>`` anywhere
    in it -- the resolver owns every one of them.
    """

    import ast
    from pathlib import Path

    import chemsmart.settings.submitters as submitters_module

    source_path = Path(submitters_module.__file__)
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    owned = {
        "num_cores",
        "mem_gb",
        "num_gpus",
        "num_nodes",
        "num_threads",
    }
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Attribute) or node.attr not in owned:
            continue
        value = node.value
        if (
            isinstance(value, ast.Attribute)
            and value.attr == "server"
            and isinstance(value.value, ast.Name)
            and value.value.id == "self"
        ):
            offenders.append(f"{source_path.name}:{node.lineno} .{node.attr}")

    assert not offenders, (
        "these read a resource straight off the server profile instead of "
        "the resolved request, so the ceiling does not bound them: "
        f"{offenders}"
    )


def test_the_envelope_s_episode_window_is_the_wall_clock_asked_for(tmp_path):
    """The job runs one whole cycle, so the bound is the episode window
    plus the postprocessing reserve the human granted."""

    from types import SimpleNamespace as NS

    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = _slurm_server()  # NUM_HOURS=24
    envelope = NS(
        episode_wall_time_seconds=5400.0,  # the live sm3-run3 envelope
        postprocess_reserve_seconds=600.0,
    )
    request = resolve_scheduler_request(
        resources=_resources(cores=4, memory_gb=8),
        server=server,
        sealed=True,
        envelope=envelope,
    )
    # 6000 s -> 2 h, not the profile's 24.
    assert request.requested_hours == 2
    assert request.hours == 2
    assert not request.clamped


def test_the_profile_caps_a_wall_clock_it_will_not_grant(tmp_path):
    """Wall time is capped, not refused: a shorter clock does not change
    the method, and the host's own node timeout sits far below both.
    The asymmetry with cores and memory is deliberate."""

    from types import SimpleNamespace as NS

    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = _slurm_server()  # NUM_HOURS=24
    envelope = NS(
        episode_wall_time_seconds=400000.0,  # ~111 h
        postprocess_reserve_seconds=0.0,
    )
    request = resolve_scheduler_request(
        resources=_resources(cores=4, memory_gb=8),
        server=server,
        sealed=True,
        envelope=envelope,
    )
    assert request.requested_hours == 112
    assert request.hours == 24
    assert request.clamped
    assert any("wall time capped" in line for line in request.observations)


def test_an_unsealed_run_is_the_escape_from_a_profile_too_small(tmp_path):
    """A 4 GB profile must not be a dead end for a 2 GB calculation.

    The sealed ceiling subtracts the headroom before it looks at the
    request, so on any profile declaring <= 6 GB -- a small VM, a CI
    runner, a container queue -- every sealed request was refused,
    including ones the machine could trivially satisfy. --unsealed drops
    the headroom and holds the request to the profile alone.
    """

    from chemsmart.agent._contracts import ContractError
    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    small = Server(
        "small",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=2,
        MEM_GB=4,
        NUM_HOURS=1,
        QUEUE_NAME="compute",
    )
    request = _resources(cores=1, memory_gb=2)

    with pytest.raises(ContractError, match="headroom"):
        resolve_scheduler_request(resources=request, server=small, sealed=True)

    unsealed = resolve_scheduler_request(
        resources=request, server=small, sealed=False
    )
    assert (unsealed.cores, unsealed.memory_gb) == (1, 2)
    assert not unsealed.clamped


def test_the_sealed_flag_reaches_the_dispatcher_from_the_command_line():
    """A flag nothing threads is a flag that does not exist."""

    import inspect

    from chemsmart.agent.driver import GoalDriver

    assert "sealed" in inspect.signature(GoalDriver.__init__).parameters
    assert (
        inspect.signature(GoalDriver.__init__).parameters["sealed"].default
        is True
    )

    from chemsmart.cli.agent import goal

    flags = {
        option
        for parameter in goal.params
        for option in getattr(parameter, "opts", ())
    }
    assert "--unsealed" in flags or "--sealed" in flags
    # And it is persisted, so a resumed goal keeps the ceiling it ran
    # under rather than silently changing it on the next cycle.
    from chemsmart.agent.driver import GoalDriver as _D

    source = inspect.getsource(_D.__init__)
    assert '"sealed": bool(sealed)' in source
