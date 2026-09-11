"""A plan with no calculation node, made at wake, is executed by the host.

A woken session held every number the goal asked for, planned a correct
seventeen-node analysis-only chain naming the declared ids, and stopped
to wait for a review nothing would ever build; the goal returned to the
human with sixteen engine calls unspent (NOVEL-2 po2, 2026-09-04). The
owner ruled such a plan is admitted and executed with no displayed
review: it changes no identity, state or condition and launches no
engine. These pins say the executor's analysis phase walks it without a
bundle, the session host walks it when planned under a goal wake and not
otherwise, and the driver records the revision it is.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.executor import execute_analysis_only_toolchain
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    RegisteredResultInputIntentV1,
    build_scientific_toolchain_plan,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

from .test_the_executor_walks_the_approved_analysis_chain import (
    _RESULT,
    _analysis_node,
)
from .test_the_goal_loop_recovers_or_returns import (
    _delivery_rows,
    _loop,
    _planning_session,
)

pytestmark = pytest.mark.capability("tool:plan_scientific_workflow")


def _registered_chain(registered_id):
    extraction = _analysis_node(
        "extract-registered",
        "result_extraction",
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="raw", artifact_id=registered_id
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("extract-registered",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy",
                source_kind="analysis_output",
                producer_node_id="extract-registered",
                producer_output_id="e",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy",
                quantity_kind="energy",
                unit="hartree",
            ),
        ),
    )
    return build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(),
        calculation_observables={},
        analysis_nodes=(extraction, claims),
        required_output_ids=("final-energy",),
    )


def _host_with_result(tmp_path, **extra):
    resolved = _RESULT.resolve()
    artifact = TrustedArtifactRefV1(
        artifact_id="orca-result-8ae1cdc683f8eb7d",
        kind="orca_output",
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "run" / "events.jsonl", session_id="wake"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
        **extra,
    )
    host.artifacts[artifact.artifact_id] = artifact
    return host, artifact.artifact_id


def test_the_analysis_phase_walks_a_chain_with_no_bundle(tmp_path):
    host, registered_id = _host_with_result(tmp_path)
    toolchain = _registered_chain(registered_id)
    record = execute_analysis_only_toolchain(
        host=host,
        toolchain=toolchain,
        run_directory=tmp_path / "run",
        task_spec_sha256="a" * 64,
        workspace=tmp_path / "workspace",
    )
    assert record["analysis_status"] == "completed"
    assert record["engine_calls_consumed"] == 0
    assert {node["state"] for node in record["executed_nodes"]} == {"executed"}
    assert len(record["completion_receipt_sha256s"]) == 1
    assert "final-energy" in Path(record["report_path"]).read_text()
    kinds = [event.kind for event in host.event_store.read_events()]
    assert "analysis_claims_recorded" in kinds
    assert "analysis_completion_evaluated" in kinds
    assert kinds[-1] == "analysis_only_plan_executed"


def test_the_host_walks_an_analysis_only_plan_only_under_a_wake(tmp_path):
    host, registered_id = _host_with_result(
        tmp_path,
        execute_analysis_only_plans=True,
        analysis_only_run_directory=tmp_path / "run",
        analysis_only_workspace=tmp_path / "workspace",
    )
    toolchain = _registered_chain(registered_id)
    record = host._execute_analysis_only_plan(toolchain)
    assert record["analysis_status"] == "completed"
    assert "standing decision" in record["meaning"]

    quiet, registered_id = _host_with_result(tmp_path / "quiet")
    assert quiet.execute_analysis_only_plans is False


def test_a_chain_naming_a_calculation_is_not_walked_this_way(tmp_path):
    from chemsmart.agent._contracts import ContractError

    host, _ = _host_with_result(tmp_path)
    with pytest.raises(ContractError, match="no calculation node"):
        execute_analysis_only_toolchain(
            host=host,
            toolchain=SimpleNamespace(
                calculation_node_ids=("sp",),
                analysis_nodes=(),
                workflow_id="w",
                plan_sha256="b" * 64,
            ),
            run_directory=tmp_path / "run",
            task_spec_sha256="a" * 64,
            workspace=tmp_path,
        )


def test_the_driver_records_the_walk_as_the_revision_it_is(tmp_path):
    rows = list(_delivery_rows()) + [
        {
            "kind": "analysis_only_plan_executed",
            "payload": {
                "toolchain_plan_sha256": "e" * 64,
                "analysis_status": "completed",
                "engine_calls_consumed": 0,
            },
        }
    ]
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", terminal="planned", wake_rows=rows),
        ],
        executes=[],
    )
    assert result.settlement == "achieved"
    ledger_path = (
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    ) / "ledger.jsonl"
    entries = [
        json.loads(line)
        for line in ledger_path.read_text(encoding="utf-8").splitlines()
    ]
    admitted = [e for e in entries if e["kind"] == "revision_admitted"]
    assert len(admitted) == 1
    assert admitted[0]["payload"]["analysis_only"] is True
    assert admitted[0]["payload"]["checks"]["engine_calls_within_budget"]
    assert result.revisions_admitted == 1


def _chain_with_a_planned_estimator(registered_id):
    """The same chain, plus the expression that computes the spread."""

    extraction = _analysis_node(
        "extract-registered",
        "result_extraction",
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="raw", artifact_id=registered_id
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    # A second reading of the same result, so the estimator spans two
    # receipts. That was required when a single-receipt spread was
    # refused; it is now reported rather than refused, and the chain is
    # kept this shape because it is what a real spread looks like.
    second = _analysis_node(
        "extract-again",
        "result_extraction",
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="raw2", artifact_id=registered_id
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e2", selector="scf_energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e2", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    spread = _analysis_node(
        "spread",
        "quantity_expression",
        dependencies=("extract-again", "extract-registered"),
        inputs=(
            AnalysisInputIntentV1(
                input_id="e2_in",
                source_kind="analysis_output",
                producer_node_id="extract-again",
                producer_output_id="e2",
            ),
            AnalysisInputIntentV1(
                input_id="e_in",
                source_kind="analysis_output",
                producer_node_id="extract-registered",
                producer_output_id="e",
            ),
        ),
        expression_nodes=(
            # Host vocabulary over a computed value and no coefficient
            # of the session's own: an estimator built from a
            # model-authored scale factor or exponent is the model's
            # number, and the provenance walk refuses it as evidence.
            # What survives is what a spread actually is -- a
            # difference, a max minus a min, an absolute value over
            # numbers the host computed.
            {
                "node_id": "diff",
                "operation": "subtract",
                "input_ids": ("e_in", "e2_in"),
            },
            {
                "node_id": "half",
                "operation": "abs",
                "input_ids": ("diff",),
            },
        ),
        expression_output_node_ids=("half",),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="half", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("extract-registered", "spread"),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy",
                source_kind="analysis_output",
                producer_node_id="extract-registered",
                producer_output_id="e",
                uncertainty_producer_node_id="spread",
                uncertainty_producer_output_id="half",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy",
                quantity_kind="energy",
                unit="hartree",
            ),
        ),
    )
    return build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(),
        calculation_observables={},
        analysis_nodes=(extraction, second, spread, claims),
        required_output_ids=("final-energy",),
    )


def test_a_planned_estimator_is_evaluated_in_the_same_walk(tmp_path):
    """An estimate cannot predate its results; its estimator can.

    A planned claim node was written while the calculation was still a
    plan, so it carried no uncertainty and an executed delivery arrived
    `unstated` by construction -- costing a further cycle to say what
    the run had already computed. A plan may now name the analysis
    output that computes the spread, and the host evaluates it in the
    same provider-free walk and copies it into the claim. The executor
    names the number and never writes it, exactly as for the claimed
    value, so the claim is `measured` and cites the output it came from
    (owner ruling, 2026-09-09).
    """

    host, registered_id = _host_with_result(tmp_path)
    record = execute_analysis_only_toolchain(
        host=host,
        toolchain=_chain_with_a_planned_estimator(registered_id),
        run_directory=tmp_path / "run",
        task_spec_sha256="a" * 64,
        workspace=tmp_path / "workspace",
    )
    assert record["analysis_status"] == "completed"
    assert record["engine_calls_consumed"] == 0

    (claim_record,) = tuple(host.analysis_claim_records.values())
    (claim,) = claim_record.claims
    assert claim.uncertainty_basis == "measured"
    assert ":half" in claim.uncertainty_reference
    # The host's own number, copied from its own receipt, and non-zero
    # because it spans two readings that differ.
    assert claim.uncertainty > 0.0
