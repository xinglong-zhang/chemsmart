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


# --- The server profile is the resource authority -------------------
#
# The operator's server.yaml / SLURM.yaml decides how much machine a job
# gets, and the agent does not argue with it: one set of numbers reaches
# both the #SBATCH line and the engine, nothing is clamped, and nothing
# is refused after a goal is under way (owner ruling, 2026-09-16).
#
# Two earlier designs of mine are recorded in scheduler_request.py: the
# envelope as the request with the profile as a ceiling, and then a
# refusal above that ceiling. A refusal mid-goal is an agent-tool failure
# that teaches a session to carry workarounds for something the host
# should simply have decided.


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


def _script_for(server, request, tmp_path):
    job = SimpleNamespace(label="goal-g1-cycle-1", PROGRAM=None)
    return build_dispatch_script(
        submitter=server.get_submitter(job, scheduler_request=request),
        python="/opt/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )


def test_the_profile_decides_the_allocation_whatever_the_episode_asked(
    tmp_path,
):
    """The envelope is recorded, not obeyed."""

    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = _slurm_server()  # NUM_CORES=6, MEM_GB=52, NUM_HOURS=24
    request = resolve_scheduler_request(
        resources=_resources(cores=4, memory_gb=8), server=server, sealed=True
    )
    assert (request.cores, request.memory_gb) == (6, 46)
    assert request.hours == 24
    # What the episode asked for rides along as an observation.
    assert request.envelope_cores == 4
    assert request.envelope_memory_gb == 8
    assert any("approved episode asked for" in o for o in request.observations)

    lines = _sbatch_lines(_script_for(server, request, tmp_path))
    assert "#SBATCH --nodes=1 --ntasks-per-node=6 --mem=46G" in lines, lines
    assert "#SBATCH --partition=compute" in lines
    assert "#SBATCH --time=24:00:00" in lines


def test_an_episode_larger_than_the_machine_is_not_refused(tmp_path):
    """A goal already under way is never failed for asking too much.

    Submitting and then refusing is an agent-tool failure: it teaches a
    session to build control logic for a decision the host owns.
    """

    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = _slurm_server()
    request = resolve_scheduler_request(
        resources=_resources(cores=512, memory_gb=4096),
        server=server,
        sealed=True,
    )
    assert (request.cores, request.memory_gb) == (6, 46)
    lines = _sbatch_lines(_script_for(server, request, tmp_path))
    assert "#SBATCH --nodes=1 --ntasks-per-node=6 --mem=46G" in lines


def test_a_tiny_profile_is_never_a_dead_end(tmp_path):
    """No profile size refuses a dispatch.

    A 4 GB profile leaves nothing above the 6 GB headroom. That is a
    reason to take the profile's memory whole, not a reason to refuse
    every job on the host.
    """

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
    request = resolve_scheduler_request(
        resources=_resources(cores=1, memory_gb=2), server=small, sealed=True
    )
    assert (request.cores, request.memory_gb) == (2, 4)
    assert request.memory_headroom_gb == 0
    lines = _sbatch_lines(_script_for(small, request, tmp_path))
    assert "#SBATCH --nodes=1 --ntasks-per-node=2 --mem=4G" in lines


def test_the_sealed_headroom_is_the_only_thing_sealed_changes(tmp_path):
    from chemsmart.settings.scheduler_request import (
        SEALED_MEMORY_HEADROOM_GB,
        resolve_scheduler_request,
    )

    server = _slurm_server()
    sealed = resolve_scheduler_request(
        resources=None, server=server, sealed=True
    )
    unsealed = resolve_scheduler_request(
        resources=None, server=server, sealed=False
    )
    assert unsealed.memory_gb == 52
    assert sealed.memory_gb == 52 - SEALED_MEMORY_HEADROOM_GB
    assert (sealed.cores, sealed.hours) == (unsealed.cores, unsealed.hours)


def test_a_fractional_memory_is_named_in_a_unit_that_can_hold_it(tmp_path):
    """`#SBATCH --mem=87.5G` is not a specification Slurm accepts.

    The live CUHK profile declares MEM_GB 160 while its own wizard comment
    records 95777 MB (93.5 GB) per node. An operator correcting that makes
    the sealed memory 87.5 and every dispatch on the host unsubmittable.
    Only the integer was masking it.
    """

    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    server = Server(
        "cuhk-true-memory",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=64,
        MEM_GB=93.5,
        NUM_HOURS=24,
        QUEUE_NAME="chpc",
    )
    request = resolve_scheduler_request(
        resources=None, server=server, sealed=True
    )
    assert request.memory_gb == pytest.approx(87.5)
    assert request.memory_directive == "89600M"
    assert 89600 / 1024 == pytest.approx(87.5)
    line = next(
        line
        for line in _sbatch_lines(_script_for(server, request, tmp_path))
        if "--mem=" in line
    )
    assert "--mem=89600M" in line, line


def test_whole_gigabytes_keep_the_spelling_they_always_had():
    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    request = resolve_scheduler_request(
        resources=None, server=_slurm_server(), sealed=False
    )
    assert request.memory_directive == "52G"


def test_a_submitter_with_no_request_resolves_the_profile_unchanged(tmp_path):
    """`chemsmart sub` is byte-identical: no envelope, no headroom."""

    server = _slurm_server()
    job = SimpleNamespace(label="j", PROGRAM=None)
    from chemsmart.settings.scheduler_request import resolve_scheduler_request

    explicit = build_dispatch_script(
        submitter=server.get_submitter(
            job,
            scheduler_request=resolve_scheduler_request(
                resources=None, server=server, sealed=False
            ),
        ),
        python="/p",
        approval_file=tmp_path / "b.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    implicit = build_dispatch_script(
        submitter=server.get_submitter(job),
        python="/p",
        approval_file=tmp_path / "b.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        goal_id="g1",
    )
    assert explicit == implicit


def test_the_resolver_is_the_only_author_of_a_scheduler_resource_line():
    """No submitter may read a resource off the profile directly.

    Not because the profile is wrong -- it is the authority -- but because
    one resolved object is where the headroom and the fractional-memory
    spelling are applied, and a writer that goes around it gets neither.
    """

    import ast
    from pathlib import Path

    import chemsmart.settings.submitters as submitters_module

    source_path = Path(submitters_module.__file__)
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    owned = {"num_cores", "mem_gb", "num_gpus", "num_nodes", "num_threads"}
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
        f"the resolved request: {offenders}"
    )


def test_the_allocation_reaches_the_goal_s_own_record(tmp_path, monkeypatch):
    """The ledger is what a later process reads."""

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
    assert request is not None
    assert request["applied"] == {
        "cores": 6,
        "memory_gb": 46,
        "gpu_count": 0,
        "hours": 24,
    }
    assert request["source"] == "server_profile"
    assert request["envelope"] == {"cores": 4, "memory_gb": 8}

    from chemsmart.agent.driver import _dispatch_ledger_keys, _record_of

    carried = {
        key: value
        for key, value in _record_of(receipt).items()
        if key in _dispatch_ledger_keys()
    }
    assert carried["scheduler_request"]["applied"]["cores"] == 6
    assert carried["wake_command"].endswith("--goal g1")


def test_the_sealed_flag_reaches_the_dispatcher_from_the_command_line():
    """A flag nothing threads is a flag that does not exist."""

    import inspect

    from chemsmart.agent.driver import GoalDriver

    parameters = inspect.signature(GoalDriver.__init__).parameters
    assert "sealed" in parameters and parameters["sealed"].default is True

    from chemsmart.cli.agent import goal

    flags = {
        option
        for parameter in goal.params
        for option in getattr(parameter, "opts", ())
    }
    assert "--unsealed" in flags or "--sealed" in flags
    assert '"sealed": bool(sealed)' in inspect.getsource(GoalDriver.__init__)
