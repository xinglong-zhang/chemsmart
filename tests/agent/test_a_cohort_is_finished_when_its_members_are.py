"""'Finished' is a question about the wave, not about a file existing.

`wait_for_dispatched_run` and `agent wake` both asked
`(run_directory / "execution-result.json").is_file()`. The job script
creates that file with a shell redirect **before** `agent run` starts, so
it is a file from second zero -- for three hours of a running
calculation, and forever after a job killed before it wrote anything.
That is already wrong at N=1.

Under a wave it is worse: N elements redirect to one path, so the last
writer defines the cycle, each element reports `partial` by construction
(no element ever walks every approved node), and a shorter record written
after a longer one leaves the tail of the first behind it.

So each element writes its own result, and the barrier asks the durable
stream: has every member of this cohort reached a terminal state?
"""

from __future__ import annotations

import json

from chemsmart.agent.cohort import (
    build_cohort_manifest,
    cohort_completion,
    execution_result_file,
)


def _manifest(tmp_path, nodes=("a1", "a2", "a3")):
    manifest = build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="a" * 64,
        node_ids=nodes,
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    )
    manifest.write(tmp_path)
    return manifest


def _events(tmp_path, rows):
    (tmp_path / "events.jsonl").write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )


def _observed(node_id):
    return {
        "kind": "program_execution_observed",
        "payload": {"node_id": node_id, "record": {"node_id": node_id}},
    }


def _reserved(node_id, reserved_at, lease_seconds):
    return {
        "kind": "workflow_node_launch_reserved",
        "payload": {
            "node_id": node_id,
            "record": {
                "node_id": node_id,
                "reserved_at": reserved_at,
                "lease_seconds": lease_seconds,
            },
        },
    }


def test_each_element_writes_its_own_result(tmp_path):
    assert execution_result_file(tmp_path, element=None).name == (
        "execution-result.json"
    )
    assert execution_result_file(tmp_path, element=0).name == (
        "execution-result.0.json"
    )
    assert execution_result_file(tmp_path, element=2).name == (
        "execution-result.2.json"
    )


def test_a_wave_is_unfinished_while_a_member_has_not_ended(tmp_path):
    _manifest(tmp_path)
    _events(tmp_path, [_observed("a1"), _observed("a2")])
    complete, pending = cohort_completion(tmp_path)
    assert not complete
    assert pending == ("a3",)


def test_a_wave_is_finished_when_every_member_has_ended(tmp_path):
    _manifest(tmp_path)
    _events(tmp_path, [_observed("a1"), _observed("a2"), _observed("a3")])
    complete, pending = cohort_completion(tmp_path)
    assert complete
    assert pending == ()


def test_a_member_inside_its_lease_keeps_the_wave_open(tmp_path):
    """A reservation is not a result. An element still running holds the
    barrier even though its node has been touched."""

    from datetime import datetime, timedelta, timezone

    fresh = (datetime.now(timezone.utc) - timedelta(seconds=5)).isoformat()
    _manifest(tmp_path)
    _events(
        tmp_path,
        [
            _observed("a1"),
            _observed("a2"),
            _reserved("a3", fresh, 900),
        ],
    )
    complete, pending = cohort_completion(tmp_path)
    assert not complete
    assert pending == ("a3",)


def test_an_empty_result_file_does_not_mean_finished(tmp_path):
    """The shell truncates the redirect target before the engine starts.

    `is_file()` was true from second zero, so a wake fired over a run
    that had not begun.
    """

    _manifest(tmp_path)
    _events(tmp_path, [_observed("a1")])
    (tmp_path / "execution-result.0.json").write_text("", encoding="utf-8")
    (tmp_path / "execution-result.json").write_text("", encoding="utf-8")
    complete, pending = cohort_completion(tmp_path)
    assert not complete
    assert set(pending) == {"a2", "a3"}


def test_without_a_cohort_the_question_is_not_asked(tmp_path):
    """A single-job dispatch keeps the behaviour it had."""

    complete, pending = cohort_completion(tmp_path)
    assert complete is None
    assert pending == ()


def _state_changed(node_id, state):
    return {
        "kind": "workflow_node_state_changed",
        "payload": {
            "node_id": node_id,
            "node_state": state,
            # The run's own record, whose "state" is the
            # run summary and never this node's word. It is
            # deliberately different, so a reader that takes
            # it goes red instead of agreeing by accident.
            "record": {"node_id": node_id, "state": "running"},
        },
    }


def test_a_cancelled_member_ends_the_wave_rather_than_holding_it(tmp_path):
    """The barrier is terminality, not an execution receipt.

    The first predicate asked "does every member have a
    program_execution_observed event". A member cancelled before launch,
    or refused admission, reaches a terminal state without ever running
    an engine -- so the wave waited on it forever, and a cohort
    containing one cancelled calculation could never wake the Agent.
    """

    _manifest(tmp_path)
    _events(
        tmp_path,
        [
            _observed("a1"),
            _observed("a2"),
            _state_changed("a3", "cancelled"),
        ],
    )
    complete, pending = cohort_completion(tmp_path)
    assert complete, f"a cancelled member held the wave open: {pending}"


def test_a_member_still_pending_is_not_mistaken_for_terminal(tmp_path):
    _manifest(tmp_path)
    _events(tmp_path, [_observed("a1"), _state_changed("a2", "running")])
    complete, pending = cohort_completion(tmp_path)
    assert not complete
    assert set(pending) == {"a2", "a3"}


def test_a_damaged_manifest_does_not_silently_remove_the_wave(tmp_path):
    """`None` means "no cohort", and a truncated file must not say that.

    `cohort_frontier(ready, None)` admits every ready node, so a manifest
    that failed to parse would silently turn a bounded wave back into the
    flowing walk it exists to prevent -- executing work the Agent did not
    ask for in this wave, with nothing on disk saying why.
    """

    import pytest

    from chemsmart.agent._contracts import ContractError
    from chemsmart.agent.cohort import (
        COHORT_MANIFEST_FILE,
        read_cohort_manifest,
    )

    _manifest(tmp_path)
    (tmp_path / COHORT_MANIFEST_FILE).write_text(
        '{"schema_version": "chemsmart.coho', encoding="utf-8"
    )
    with pytest.raises(ContractError, match="unreadable|damaged"):
        read_cohort_manifest(tmp_path)

    (tmp_path / COHORT_MANIFEST_FILE).write_text("[]", encoding="utf-8")
    with pytest.raises(ContractError, match="unreadable|damaged"):
        read_cohort_manifest(tmp_path)


def test_a_cohort_is_written_once_and_not_replaced(tmp_path):
    """Membership is fixed at dispatch.

    `write()` overwrote unconditionally, so a second correctly-digested
    manifest replaced the first through the public writer -- a different
    experiment than the one the barrier is waiting for.
    """

    import pytest

    from chemsmart.agent._contracts import ContractError

    _manifest(tmp_path, nodes=("a1", "a2", "a3"))
    with pytest.raises(ContractError, match="already"):
        _manifest(tmp_path, nodes=("b1", "b2"))
