"""The wave contract holds on the default dispatch too.

The stem tells every session, without qualification: *"The host runs them
concurrently and wakes you once, when every member has reached a terminal
state... A node whose dependency clears mid-wave is not started for you."*
`select_execution_wave`'s own reply says *"this wave is what will be
submitted."*

On `--dispatch local` -- the default, and what the terminal interface
uses -- no manifest was written anywhere, so `_cohort_scope()` was
`None`, the walk was unbounded, and the executor ran the whole approved
DAG straight past the barrier. Nothing told the model which dispatch mode
it was in, so the rule was true of one path and asserted on both, which
is this round's own defect class in the prose rather than the code.

A local run is one process rather than N, so concurrency is genuinely not
on offer -- and concurrency was never the promise. The promise is that a
wave is what runs and then the session reasons again, and one process can
keep that exactly.
"""

from __future__ import annotations

import json
from types import SimpleNamespace


from chemsmart.agent.cohort import read_cohort_manifest


def _driver(tmp_path, *, wave, execute):
    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import (
        _envelope_file,
        _planning_session,
        _review_payload,
    )

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    (tmp_path / "bundle.json").write_text(
        json.dumps({"bundle_sha256": "7" * 64}), encoding="utf-8"
    )

    def plan_session(**kw):
        session = _planning_session("live-1", review=_review_payload())(
            workspace, kw
        )
        object.__setattr__(session, "selected_execution_wave", tuple(wave))
        return session

    return GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="g1",
        granted_by="tester",
        plan_session=plan_session,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=execute,
    )


def test_a_local_run_writes_the_manifest_the_walk_reads(tmp_path):
    """One process, still bounded to the wave the Agent chose."""

    seen: dict = {}

    def execute(*, approval_file, workspace, run_directory, **_kw):
        seen["run_directory"] = run_directory
        seen["manifest"] = read_cohort_manifest(run_directory)
        return SimpleNamespace(status="completed", analysis_status="")

    driver = _driver(tmp_path, wave=("opt-a", "opt-b"), execute=execute)
    driver.run()

    manifest = seen.get("manifest")
    assert manifest is not None, (
        "a local run of a selected wave wrote no manifest, so the "
        "executor walked the whole approved DAG past the barrier the "
        "stem promises every session"
    )
    assert tuple(manifest.node_ids) == ("opt-a", "opt-b")
    # Written before the executor is handed the directory, for the same
    # reason the scheduler path writes it before submitting.
    assert manifest.cycle == 1


def test_a_local_run_without_a_selected_wave_is_unchanged(tmp_path):
    """Every goal that predates waves keeps its own walk."""

    seen: dict = {}

    def execute(*, approval_file, workspace, run_directory, **_kw):
        seen["manifest"] = read_cohort_manifest(run_directory)
        return SimpleNamespace(status="completed", analysis_status="")

    driver = _driver(tmp_path, wave=(), execute=execute)
    driver.run()
    assert seen.get("manifest") is None
