"""The engine budget counts the engine; the host's reading is timed beside it.

REACH-1 po3 (2026-09-06): the executor stamped a node's finish after
its own evaluation, so 74 minutes of parser time on one saddle search
was charged to the goal's engine-wall grant (receipt window 15 267 s
against ORCA's 10 801 s), and the binding budget line the next session
quoted was inflated by 41 %.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from chemsmart.agent.execution import build_program_execution_receipt
from tests.agent.test_a_terminal_state_is_derived_not_grepped import (  # noqa: F401
    _reserve,
)


def _receipt(invocation, **timing):
    return build_program_execution_receipt(
        invocation,
        execution_state="engine_complete",
        exit_status=0,
        child_exit_status=0,
        engine_complete=True,
        validated=False,
        started_at="2026-09-06T12:00:00+00:00",
        finished_at="2026-09-06T12:00:10+00:00",
        **timing,
    )


@pytest.mark.capability("rule:terminal_state_vocabulary")
def test_an_old_receipt_keeps_its_digest_and_a_new_one_carries_host_time(
    tmp_path,
):
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.terminal_states import (
        derive_run_outcome,
        read_run_events,
    )

    store = RuntimeEventStore(tmp_path / "events.jsonl", session_id="s")
    _, plan, _m, _a, invocation = _reserve(store, tmp_path)
    plain = _receipt(invocation)
    assert plain.evaluated_at == ""
    # The digest is computed without the empty field, exactly as before.
    assert replace(plain, receipt_sha256=plain.receipt_sha256)
    timed = _receipt(invocation, evaluated_at="2026-09-06T12:05:10+00:00")
    assert timed.receipt_sha256 != plain.receipt_sha256
    store.record_program_execution_receipt(
        turn_id="turn-1",
        workflow_id=plan.workflow_id,
        run_id="run.water-approval",
        receipt=timed,
    )
    outcome = derive_run_outcome(read_run_events(tmp_path / "events.jsonl"))
    (node,) = outcome.nodes
    assert node.wall_seconds == pytest.approx(10.0)
    assert node.host_seconds == pytest.approx(300.0)
    assert outcome.engine_wall_seconds == pytest.approx(10.0)
    assert outcome.host_seconds == pytest.approx(300.0)
    assert outcome.public_record()["host_seconds"] == pytest.approx(300.0)
