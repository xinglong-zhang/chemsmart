"""A parked cycle's workspace record must carry the level it ran at.

`record_run` reads the level of every result from the cycle's displayed
review, through `self.review_file`. That attribute is set in `_plan` --
and a parked cycle's workspace record is written by the *wake*, which
`resume` rebuilds at the outcome phase without planning. So
`review_file` was `None` for every scheduler-dispatched cycle, and every
result row went to disk with no level and every claim with no ancestry
to resolve.

Measured on the live CUHK run (goal `butane-wave-2`): the review carries
`project_settings_text_sha256` and the full settings text for all six
nodes, and the workspace record's seven result rows all carry an empty
`level_sha256` -- so the per-claim level attribution this round built had
nothing to attribute. The machinery was right and the data never reached
it, which is the defect class in its quietest form: nothing failed.

The path is deterministic and the driver already computes it in `_plan`.
"""

from __future__ import annotations

import json



def _goal_with_review(tmp_path):
    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="g1",
        granted_by="tester",
    )
    reviews = driver.goal_dir / "reviews"
    reviews.mkdir(parents=True, exist_ok=True)
    (reviews / "cycle-2.json").write_text(
        json.dumps(
            {
                "workflow_execution_review": {
                    "node_reviews": [
                        {
                            "node_id": "opt-anti",
                            "project_settings_text": '{"basis":"6-31g*"}',
                            "project_settings_text_sha256": "a" * 64,
                        }
                    ]
                }
            }
        ),
        encoding="utf-8",
    )
    return driver


def test_a_resumed_cycle_finds_its_own_review(tmp_path):
    driver = _goal_with_review(tmp_path)
    driver.cycles = 2
    assert driver.review_file is None

    resolved = driver._review_file_for_cycle()
    assert resolved is not None
    assert resolved.name == "cycle-2.json"
    assert resolved.is_file()


def test_the_workspace_record_is_given_that_review(tmp_path):
    """Drive the real recorder: the level must reach the row."""

    from chemsmart.agent.workspace_record import read_workspace_record

    driver = _goal_with_review(tmp_path)
    driver.cycles = 2
    driver.run_directory = driver.goal_dir / "runs" / "cycle-2"
    driver.run_directory.mkdir(parents=True)
    events = driver.run_directory / "events.jsonl"
    events.write_text(
        json.dumps(
            {
                "kind": "program_result_verified",
                "payload": {
                    "record": {
                        "node_id": "opt-anti",
                        "state": "valid",
                        "jobtype": "opt",
                        "observations": {"jobtype": "opt", "pyscf": {}},
                        "output_artifacts": [{"sha256": "c" * 64}],
                    }
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    driver._record_workspace(events, "goals/g1/runs/cycle-2")

    rows = [
        row
        for row in read_workspace_record(driver.workspace)
        if row.get("kind") == "result"
    ]
    assert rows, "nothing was projected"
    assert rows[0]["level_sha256"] == "a" * 64, (
        "the woken cycle recorded its result with no level, so every "
        "claim standing on it has no ancestry to resolve"
    )


def test_a_cycle_with_no_review_on_disk_records_as_it_did(tmp_path):
    """Absence stays absence; nothing is invented."""

    driver = _goal_with_review(tmp_path)
    driver.cycles = 9
    assert driver._review_file_for_cycle() is None
