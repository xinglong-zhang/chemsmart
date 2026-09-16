"""The wake asks the barrier, not the filesystem.

`cohort_completion` is the host's own answer to "has this wave ended" --
terminality derived from the durable stream, a member still inside its
launch lease excluded, a cancelled member counted as ended. It was
written, tested against itself, and **read by nothing**: the declared-
but-no-reader defect this round exists to remove, in the one place where
it decides whether the Agent is ever woken.

What actually gated the wake was `(run_directory / "execution-result.json")
.is_file()` -- the predicate §2 row 3 of the plan falsified, because the
job script's own `>` redirect creates that file before the engine starts.
For a cohort it is worse than early: the file is never created at all,
because every element writes `execution-result.<i>.json`. So the wake job
fires, `chemsmart agent wake` refuses with "has not written its result
yet", and the goal parks forever with all its science finished on disk.

Two authors for one path is the other half: `execution_result_file` owns
the per-element name and the array script spelled it inline.
"""

from __future__ import annotations

import json

from click.testing import CliRunner

from chemsmart.agent.cohort import (
    build_cohort_manifest,
    cohort_completion,
    execution_result_file,
)


def _manifest(run_directory, node_ids=("a1", "a2")):
    run_directory.mkdir(parents=True, exist_ok=True)
    build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="e" * 64,
        node_ids=node_ids,
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    ).write(run_directory)


def _stream(run_directory, terminal):
    rows = [
        {
            "kind": "workflow_node_state_changed",
            "payload": {
                "node_id": node,
                "node_state": state,
                "record": {"node_id": node, "state": "running"},
            },
        }
        for node, state in terminal.items()
    ]
    (run_directory / "events.jsonl").write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )


def test_the_array_script_writes_the_path_the_host_owns(tmp_path):
    """One author for the per-element result path."""

    from types import SimpleNamespace

    from chemsmart.agent.dispatch import build_cohort_dispatch_script
    from chemsmart.settings.server import Server

    server = Server(
        "canned-slurm",
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=4,
        MEM_GB=8,
        NUM_HOURS=1,
        QUEUE_NAME="compute",
    )
    submitter = server.get_submitter(
        SimpleNamespace(label="goal-g1-cycle-1", PROGRAM=None, folder=".")
    )
    script = build_cohort_dispatch_script(
        submitter=submitter,
        python="/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "run",
        cohort_size=2,
    )
    expected = str(execution_result_file(tmp_path / "run", element=0)).replace(
        "execution-result.0.json", "execution-result."
    )
    assert expected in script, (
        "the array script spells the per-element result path itself, so "
        "the host has two authors for where an element writes"
    )


def test_a_cohort_with_a_live_member_is_not_over(tmp_path):
    run_directory = tmp_path / "run"
    _manifest(run_directory)
    _stream(run_directory, {"a1": "validated"})
    complete, pending = cohort_completion(run_directory)
    assert complete is False
    assert pending == ("a2",)


def test_a_cohort_whose_members_all_ended_is_over(tmp_path):
    run_directory = tmp_path / "run"
    _manifest(run_directory)
    _stream(run_directory, {"a1": "validated", "a2": "cancelled"})
    complete, pending = cohort_completion(run_directory)
    assert complete is True, pending


def test_the_wake_command_asks_the_barrier_for_a_cohort(tmp_path, monkeypatch):
    """Drive the public command: it must not demand the single-job file.

    A cohort's elements each write `execution-result.<i>.json`, so the
    file the wake insisted on is one nothing ever creates. Every member
    has ended and every receipt is on disk, and the goal cannot re-enter.
    """

    from chemsmart.agent import driver as driver_module
    from chemsmart.agent.dispatch import EXECUTION_RESULT_FILE
    from chemsmart.cli.agent import agent

    workspace = tmp_path / "ws"
    run_directory = workspace / ".chemsmart-agent" / "goals" / "g1" / "runs"
    run_directory = run_directory / "cycle-1"
    _manifest(run_directory, ("a1", "a2"))
    _stream(run_directory, {"a1": "validated", "a2": "failed"})
    # Every element wrote its own result; the single-job name is absent.
    for element in (0, 1):
        execution_result_file(run_directory, element=element).write_text(
            json.dumps({"status": "completed"}), encoding="utf-8"
        )
    assert not (run_directory / EXECUTION_RESULT_FILE).exists()

    ran: list[str] = []

    class _Resumed:
        cycles = 1

        def __init__(self):
            self.run_directory = run_directory

        def run(self):
            ran.append("ran")
            return driver_module.GoalLoopResultV1(
                goal_id="g1",
                settlement="achieved",
                cycles=1,
                revisions_admitted=0,
            )

    monkeypatch.setattr(
        driver_module.GoalDriver,
        "resume",
        classmethod(lambda cls, **kw: _Resumed()),
    )
    result = CliRunner().invoke(
        agent, ["wake", "--workspace", str(workspace), "--goal", "g1"]
    )
    assert ran, (
        "the wake refused a finished wave because it looked for the "
        f"single-job result file: {result.output}"
    )
    assert result.exit_code == 0, result.output


def test_the_wake_still_refuses_a_cohort_that_has_not_ended(
    tmp_path, monkeypatch
):
    """Not a weakening: an unfinished wave is still refused."""

    from chemsmart.agent import driver as driver_module
    from chemsmart.cli.agent import agent

    workspace = tmp_path / "ws"
    run_directory = (
        workspace / ".chemsmart-agent" / "goals" / "g1" / "runs" / "cycle-1"
    )
    _manifest(run_directory, ("a1", "a2"))
    _stream(run_directory, {"a1": "validated"})

    class _Resumed:
        cycles = 1

        def __init__(self):
            self.run_directory = run_directory

        def run(self):  # pragma: no cover - must not be reached
            raise AssertionError("the wake ran over an unfinished wave")

    monkeypatch.setattr(
        driver_module.GoalDriver,
        "resume",
        classmethod(lambda cls, **kw: _Resumed()),
    )
    result = CliRunner().invoke(
        agent, ["wake", "--workspace", str(workspace), "--goal", "g1"]
    )
    assert result.exit_code != 0
    assert "a2" in result.output, result.output
