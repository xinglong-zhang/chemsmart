"""How a node ended used to live in free text: a non-converged scan step
and a SHARK abort arrived as one blob, and an interrupted engine was an
empty rule-id tuple beside an English sentence nobody persisted.

The terminal vocabulary is derived on read from the sealed events plus
the artifact's own parser facts -- nothing new is written, so the
provenance of a terminal state is the event hashes and artifact digests
the derivation read. These tests drive the derivation over streams built
by the real launch fence and receipt writer, never hand-forged events.
"""

import pytest

from chemsmart.agent.execution import build_program_execution_receipt
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.terminal_states import (
    NODE_TERMINAL_STATES,
    NodeTerminalStateV1,
    derive_run_outcome,
    read_run_events,
)

from .test_runtime_v2_launch_fence import _reserve


def _stream_with_receipt(tmp_path, *, findings, execution_state="failed"):
    path = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(path, session_id="water-session")
    _, plan, _materialized, _approval, invocation = _reserve(store, tmp_path)
    receipt = build_program_execution_receipt(
        invocation,
        execution_state=execution_state,
        exit_status=1 if execution_state == "failed" else 0,
        child_exit_status=1 if execution_state == "failed" else 0,
        engine_complete=execution_state != "failed",
        validated=execution_state == "validated",
        findings=tuple(findings),
        started_at="2026-08-04T00:00:00+00:00",
        finished_at="2026-08-04T00:00:05+00:00",
    )
    store.record_program_execution_receipt(
        turn_id="turn-1",
        workflow_id=plan.workflow_id,
        run_id="run.water-approval",
        receipt=receipt,
    )
    return path


def test_a_timeout_is_its_own_ending_with_budget_facts(tmp_path):
    path = _stream_with_receipt(
        tmp_path,
        findings=(
            "execution.process.nonzero_or_unknown",
            "execution.process.timeout",
        ),
    )
    outcome = derive_run_outcome(read_run_events(path))
    assert outcome.workflow_id == "water-workflow"
    assert outcome.engine_calls_consumed == 1
    assert outcome.engine_wall_seconds == pytest.approx(5.0)
    (node,) = outcome.nodes
    assert node.state == "timeout_terminated"
    assert node.wall_seconds == pytest.approx(5.0)
    assert node.evidence_event_hashes, "the citation must not be empty"


def test_an_unconfirmed_timeout_reads_ambiguous(tmp_path):
    path = _stream_with_receipt(
        tmp_path,
        findings=(
            "execution.process.timeout",
            "execution.process.termination_ambiguous",
        ),
    )
    (node,) = derive_run_outcome(read_run_events(path)).nodes
    assert node.state == "timeout_ambiguous"


def test_a_human_interrupt_reads_as_its_own_ending(tmp_path):
    path = _stream_with_receipt(
        tmp_path, findings=("execution.process.external_signal",)
    )
    (node,) = derive_run_outcome(read_run_events(path)).nodes
    assert node.state == "external_signal_terminated"


def test_a_reservation_without_a_receipt_is_interrupted_mid_engine(
    tmp_path,
):
    path = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(path, session_id="water-session")
    _reserve(store, tmp_path)
    (node,) = derive_run_outcome(read_run_events(path)).nodes
    assert node.state == "interrupted_mid_engine"


def test_a_plain_failure_is_native_not_invented(tmp_path):
    path = _stream_with_receipt(
        tmp_path, findings=("execution.process.nonzero_or_unknown",)
    )
    (node,) = derive_run_outcome(read_run_events(path)).nodes
    assert node.state == "failed_native"
    assert node.converged is None, (
        "no artifact was readable, so convergence stays absent -- "
        "never manufactured"
    )


def test_the_scan_classification_needs_the_observed_facts():
    from chemsmart.agent.terminal_states import _classify_failure

    assert (
        _classify_failure(
            jobtype="scan",
            findings=("execution.process.nonzero_or_unknown",),
            native_class="",
            converged=False,
            reached=2,
            planned=12,
        )
        == "failed_nonconverged_scan_step"
    )
    assert (
        _classify_failure(
            jobtype="opt",
            findings=("orca.result.optimization_not_converged",),
            native_class="",
            converged=False,
            reached=None,
            planned=None,
        )
        == "failed_nonconverged_geometry"
    )
    assert (
        _classify_failure(
            jobtype="sp",
            findings=(),
            native_class="scf_convergence",
            converged=None,
            reached=None,
            planned=None,
        )
        == "failed_nonconverged_scf"
    )


@pytest.mark.capability("gate:terminal_state_vocabulary")
def test_a_recorded_cause_outranks_the_convergence_it_leaves_behind():
    """Every native class a program can report, against a dead run's flag.

    A crashed run leaves ``converged is False`` behind exactly as a
    genuinely unconverged one does, so a classification that reads the
    flag first cannot tell them apart. Only two classes are themselves
    convergence statements; the rest name a cause, and the cause wins.
    The table is walked from the program rule sets rather than from the
    two classes that were observed failing, because the next crash will
    be a third one.
    """

    from chemsmart.agent.terminal_states import (
        _CONVERGENCE_FAILURE_CLASSES,
        _UNDIAGNOSED_FAILURE_CLASSES,
        _classify_failure,
    )
    from chemsmart.io.native_failure import (
        _GAUSSIAN_RULES,
        _ORCA_RULES,
        _PYSCF_STAGE_CLASSES,
        _XTB_RULES,
    )

    classes = {
        *(name for name, _patterns in _ORCA_RULES),
        *(name for name, _patterns in _GAUSSIAN_RULES),
        *(name for name, _patterns in _XTB_RULES),
        *_PYSCF_STAGE_CLASSES.values(),
        "native_runtime",
        "incomplete_output",
        "driver_exception",
    }
    assert _CONVERGENCE_FAILURE_CLASSES < classes, (
        "the convergence classes must be drawn from the same vocabulary "
        "the programs actually report"
    )

    named = classes - _CONVERGENCE_FAILURE_CLASSES
    named -= _UNDIAGNOSED_FAILURE_CLASSES
    for native_class in sorted(named):
        for jobtype, reached, planned in (
            ("opt", None, None),
            ("scan", 2, 12),
        ):
            assert (
                _classify_failure(
                    jobtype=jobtype,
                    findings=(
                        "execution.process.nonzero_or_unknown",
                        f"orca.native_failure.{native_class}",
                        "orca.result.optimization_not_converged",
                    ),
                    native_class=native_class,
                    converged=False,
                    reached=reached,
                    planned=planned,
                )
                == "failed_native"
            ), f"{native_class} on a {jobtype} was read as a convergence "

    # And the two that are convergence statements keep their meaning, so
    # the branch above is a narrowing and not a new blanket.
    assert (
        _classify_failure(
            jobtype="opt",
            findings=("xtb.native_failure.geometry_optimization",),
            native_class="geometry_optimization",
            converged=False,
            reached=None,
            planned=None,
        )
        == "failed_nonconverged_geometry"
    )
    assert (
        _classify_failure(
            jobtype="sp",
            findings=("orca.native_failure.scf_convergence",),
            native_class="scf_convergence",
            converged=False,
            reached=None,
            planned=None,
        )
        == "failed_nonconverged_scf"
    )
    # A run that merely hit its iteration cap terminates normally and so
    # reports no class at all: the ordinary non-convergence is untouched.
    assert (
        _classify_failure(
            jobtype="opt",
            findings=("orca.result.optimization_not_converged",),
            native_class="",
            converged=False,
            reached=None,
            planned=None,
        )
        == "failed_nonconverged_geometry"
    )

    # And the classes that name no cause do not get to speak over one.
    # ORCA reports a non-converged relaxed scan step by failing to store
    # the step's geometry and aborting inside its property module, which
    # matches no rule and lands on native_runtime. Reproduced at one
    # rank with no MPI in the run, so it is neither an MPI defect nor
    # non-deterministic: the same input fails at the same step, and the
    # output says "The optimization did not converge but reached the
    # maximum number of" steps.
    for fallback in sorted(_UNDIAGNOSED_FAILURE_CLASSES):
        assert (
            _classify_failure(
                jobtype="scan",
                findings=(
                    "execution.process.nonzero_or_unknown",
                    f"orca.native_failure.{fallback}",
                    "orca.result.optimization_not_converged",
                ),
                native_class=fallback,
                converged=False,
                reached=2,
                planned=9,
            )
            == "failed_nonconverged_scan_step"
        ), fallback
        # With nothing saying anything about convergence, the same
        # fallback still ends as a native failure.
        assert (
            _classify_failure(
                jobtype="sp",
                findings=("execution.process.nonzero_or_unknown",),
                native_class=fallback,
                converged=None,
                reached=None,
                planned=None,
            )
            == "failed_native"
        ), fallback


def test_a_withdrawn_grant_survives_the_process(tmp_path):
    """The executor wrote cancelled into an in-memory record only, so
    a human's withdrawal derived as not_launched afterwards -- the
    withdrawn grant and a node that never came up shared one word.
    The durable vocabulary now carries it: a pending node cancels, a
    launched node never does, and the workflow summary word follows."""

    from chemsmart.agent.execution import (
        ContractError,
        build_workflow_run_state,
        transition_workflow_node,
    )

    from .test_runtime_v2_launch_fence import _frontier

    plan, _materialized, approval, _invocation = _frontier(tmp_path)
    run_state = build_workflow_run_state(
        plan=plan,
        approval=approval,
        run_id="run.water-approval",
        approval_consumed=True,
    )
    (row,) = run_state.nodes
    assert row.state == "pending"
    cancelled = transition_workflow_node(
        run_state,
        node_id="sp-initial",
        new_state="cancelled",
        plan=plan,
        failure_rule_ids=("execution.cancelled.human",),
        timestamp="2026-08-04T00:00:01+00:00",
    )
    (node,) = cancelled.nodes
    assert node.state == "cancelled"
    assert node.failure_rule_ids == ("execution.cancelled.human",)
    assert cancelled.state == "cancelled"
    assert cancelled.finished_at

    running = transition_workflow_node(
        run_state,
        node_id="sp-initial",
        new_state="running",
        plan=plan,
        invocation_sha256="5" * 64,
        timestamp="2026-08-04T00:00:01+00:00",
    )
    with pytest.raises(ContractError, match="invalid workflow node state"):
        transition_workflow_node(
            running,
            node_id="sp-initial",
            new_state="cancelled",
            plan=plan,
            timestamp="2026-08-04T00:00:02+00:00",
        )


def test_every_state_the_derivation_can_emit_is_declared():
    with pytest.raises(ValueError, match="unsupported node terminal"):
        NodeTerminalStateV1(
            node_id="n",
            program="orca",
            jobtype="sp",
            state="made_up",
        )
    assert "cancelled" in NODE_TERMINAL_STATES


def test_a_second_run_cannot_even_be_forged_into_a_stream(tmp_path):
    """Defense in depth, from the outside in.

    The launch fence refuses a second workflow run into a recorded
    directory outright, so a two-run stream cannot be produced by any
    production writer; and a doctored reservation event fails its own
    record digest at construction, before the run-record map could ever
    hold two entries. The derivation's exactly-one guard therefore sits
    behind two working fences -- pinned here so weakening either one
    fails a test.
    """

    import dataclasses

    from chemsmart.agent._contracts import ContractError

    path = tmp_path / "events" / "runtime.jsonl"
    store = RuntimeEventStore(path, session_id="water-session")
    _reserve(store, tmp_path)
    events = read_run_events(path)
    reservation = next(
        event
        for event in events
        if event.kind == "workflow_node_launch_reserved"
    )
    # The forge dies at event *construction*: the typed event validates
    # its reservation record's own digest, so a doctored run id cannot
    # even become an event, let alone reach the run-record map.
    with pytest.raises(ContractError, match="digest mismatch"):
        dataclasses.replace(
            reservation,
            sequence=reservation.sequence + len(events),
            payload={
                **reservation.payload,
                "run_id": "run.other",
                "record": {
                    **(reservation.payload.get("record") or {}),
                    "run_id": "run.other",
                },
            },
        )


def test_an_engine_complete_receipt_with_a_failed_validation_is_a_failure(
    tmp_path,
):
    """Observed live (R2, hydrogen peroxide from a planar start): the
    engine finished, the validator recorded result.stationary_point_order,
    the node row went to failed -- and the derivation, preferring the
    receipt's engine_complete over the row, called it unvalidated. A
    refused result is a failure with a cause, and the cause is what the
    goal driver's recovery reads."""

    from chemsmart.agent.terminal_states import STATIONARY_POINT_ORDER_FINDING

    path = _stream_with_receipt(
        tmp_path,
        findings=(STATIONARY_POINT_ORDER_FINDING,),
        execution_state="engine_complete",
    )
    (node,) = derive_run_outcome(read_run_events(path)).nodes
    assert node.state == "failed_wrong_stationary_point"

    clean = _stream_with_receipt(
        tmp_path / "clean", findings=(), execution_state="engine_complete"
    )
    (node,) = derive_run_outcome(read_run_events(clean)).nodes
    assert node.state == "engine_complete_unvalidated"


def test_an_anomaly_rides_its_node_into_the_outcome(tmp_path):
    """A host-detected surprise recorded beneath a node's verdict reaches
    the run outcome -- and so the wake context and inspect_run -- with
    its signal, status, numbers and receipt digest; a stream without
    the sensor derives an empty tuple."""

    path = _stream_with_receipt(
        tmp_path, findings=("result.stationary_point_order",)
    )
    before = derive_run_outcome(read_run_events(path))
    (node,) = before.nodes
    assert node.anomalies == ()
    store = RuntimeEventStore(path, session_id="water-session")
    store.append(
        turn_id="turn-1",
        kind="anomaly_observed",
        payload={
            "receipt_sha256": "c" * 64,
            "status": "unreplicated",
            "node_id": node.node_id,
            "signal_id": "stationary_point.unexpected_order",
            "record": {
                "signal_id": "stationary_point.unexpected_order",
                "status": "unreplicated",
                "values": {"observed_imaginary_modes": 1},
            },
        },
    )
    after = derive_run_outcome(read_run_events(path))
    (node,) = after.nodes
    (anomaly,) = node.anomalies
    assert anomaly["signal_id"] == "stationary_point.unexpected_order"
    assert anomaly["values"] == {"observed_imaginary_modes": 1}
    assert anomaly["receipt_sha256"] == "c" * 64
    assert after.public_record()["nodes"][0]["anomalies"] == (anomaly,)
