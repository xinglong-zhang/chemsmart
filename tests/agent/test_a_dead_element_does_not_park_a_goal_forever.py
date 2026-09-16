"""A member whose process died is not a member still working.

Measured live on CUHK (goal `butane-wave-1`, array 2135192): three of
four elements died before reserving a launch, so their nodes reached no
terminal state, and the barrier -- correctly, from the stream alone --
said the wave was not over. It was right and the goal was stuck: the
dependent wake job had already fired and no process would ever write
those three nodes again.

The stream cannot tell "still running" from "process gone"; that is what
the launch lease answers for a node that reserved, and what the scheduler
answers for an element that never got that far. So the wake asks the
scheduler exactly one question -- *can anything still be working on this
wave* -- and never what any calculation means. When the array has ended,
a member with no events is `not_launched`, which is already a terminal
state in the host's own vocabulary, and the Agent reads "these three
never launched" as the evidence it is.
"""

from __future__ import annotations

import json
import subprocess

from click.testing import CliRunner

from chemsmart.agent.cohort import build_cohort_manifest
from chemsmart.agent.dispatch import DISPATCH_RECEIPT_FILE, DispatchReceiptV1

_ENDED = (
    "4242_0|4242|0|COMPLETED|00:05:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:05:05\n"
    "4242_1|4242|1|FAILED|00:00:03|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:00:08\n"
)
_RUNNING = (
    "4242_0|4242|0|COMPLETED|00:05:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:05:05\n"
    "4242_1|4242|1|RUNNING|00:02:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|N/A\n"
)


def _run_directory(tmp_path):
    run_directory = (
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "g1" / "runs"
    ) / "cycle-1"
    run_directory.mkdir(parents=True)
    build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="e" * 64,
        node_ids=("a1", "a2"),
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    ).write(run_directory)
    # a1 validated; a2's element died before it reserved anything.
    (run_directory / "events.jsonl").write_text(
        json.dumps(
            {
                "kind": "workflow_node_state_changed",
                "payload": {
                    "node_id": "a1",
                    "node_state": "validated",
                    "record": {"node_id": "a1", "state": "running"},
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(
            DispatchReceiptV1(
                scheduler="SLURM",
                job_id="4242",
                submitted_at="2026-09-16T12:00:00+00:00",
                submit_command="sbatch x.sh",
                submit_script="x.sh",
                run_directory=str(run_directory),
                approval_file=str(tmp_path / "bundle.json"),
                goal_id="g1",
                cycle=1,
                wake_command="wake",
                wake_job_id="4243",
                cohort_node_ids=("a1", "a2"),
            ).public_record()
        ),
        encoding="utf-8",
    )
    return run_directory


def _wake(tmp_path, run_directory, monkeypatch, squeue_output):
    from chemsmart.agent import driver as driver_module
    from chemsmart.cli import agent as agent_cli

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
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda argv, **kw: subprocess.CompletedProcess(
            argv, 0, stdout=squeue_output
        ),
    )
    result = CliRunner().invoke(
        agent_cli.agent,
        ["wake", "--workspace", str(tmp_path / "ws"), "--goal", "g1"],
    )
    return ran, result


def test_a_wave_whose_array_has_ended_wakes_the_agent(tmp_path, monkeypatch):
    run_directory = _run_directory(tmp_path)
    ran, result = _wake(tmp_path, run_directory, monkeypatch, _ENDED)
    assert ran, (
        "the goal stayed parked over a member no process will ever write "
        f"again: {result.output}"
    )
    assert result.exit_code == 0, result.output


def test_a_wave_whose_array_is_still_running_is_still_refused(
    tmp_path, monkeypatch
):
    """Not a weakening: a live element still holds the barrier."""

    run_directory = _run_directory(tmp_path)
    ran, result = _wake(tmp_path, run_directory, monkeypatch, _RUNNING)
    assert not ran
    assert result.exit_code != 0
    assert "a2" in result.output, result.output


def test_a_wave_with_no_receipt_keeps_the_stream_s_answer(
    tmp_path, monkeypatch
):
    """Without a dispatch receipt there is no scheduler to ask."""

    run_directory = _run_directory(tmp_path)
    (run_directory / DISPATCH_RECEIPT_FILE).unlink()
    ran, result = _wake(tmp_path, run_directory, monkeypatch, _ENDED)
    assert not ran
    assert result.exit_code != 0
