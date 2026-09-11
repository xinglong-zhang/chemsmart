"""A refusal states its gate, diagnoses the cause, names a route and its
cost -- and the route is one the session can walk where it stands.

Measured over one round the agent met 470 refusals and retried the
refused tool in 60% of them: the refusal message is where the model is
taught. NOVEL-2 po2's woken session met one whose only imperative was
"plan one first", obeyed it, and lost six computed numbers.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError, RoutedContractError
from chemsmart.agent.scientific_toolchain import (
    ScientificToolchainContractError,
    build_scientific_toolchain_plan,
)

from .test_a_claim_carries_the_id_it_answers import _claim_node
from .test_an_analysis_only_plan_executes_at_wake import _host_with_result

pytestmark = pytest.mark.capability("rule:wake.refusal_is_a_deliverable")


def test_a_routed_refusal_carries_its_four_fields_and_renders_them():
    error = RoutedContractError(
        gate="g.x",
        invariant="an invariant.",
        diagnosis="what the host saw.",
        route="call a tool",
    )
    assert error.failure_report == {
        "gate": "g.x",
        "invariant": "an invariant.",
        "diagnosis": "what the host saw.",
        "route": "call a tool",
        "cost": "no engine call",
    }
    assert str(error).startswith("[g.x] an invariant. Diagnosis:")
    assert isinstance(error, ContractError)
    with pytest.raises(ContractError, match="names its gate"):
        RoutedContractError(gate="g", invariant="", diagnosis="d", route="r")


def test_an_unknown_workflow_names_the_route_where_the_session_stands(
    tmp_path,
):
    woken, _ = _host_with_result(
        tmp_path / "woken", execute_analysis_only_plans=True
    )
    with pytest.raises(RoutedContractError) as refused:
        woken._resolve_program_workflow("fccs-gauche-anti-dg")
    report = refused.value.failure_report
    assert report["gate"] == "workflow.id_names_a_plan_this_session_holds"
    assert "record_analysis_claims" in report["route"]
    assert "no calculation_nodes" in report["route"]
    assert "plan one first" not in str(refused.value)

    fresh, _ = _host_with_result(tmp_path / "fresh")
    with pytest.raises(RoutedContractError) as refused:
        fresh._resolve_program_workflow("w")
    assert "plan_scientific_workflow plans one" in (
        refused.value.failure_report["route"]
    )


def test_a_required_output_with_no_producer_names_the_route():
    with pytest.raises(ScientificToolchainContractError) as refused:
        build_scientific_toolchain_plan(
            plan_id="p",
            workflow_id="w",
            command_workflow_draft_sha256="9" * 64,
            calculation_nodes=(),
            calculation_observables={},
            analysis_nodes=(),
            required_output_ids=("cis_gs_mult",),
        )
    assert refused.value.cause == "required_output_has_no_producer"
    assert "coordinate_at_minimum" in refused.value.next_legal_route
    assert "no engine call" in refused.value.next_legal_route


def test_a_claim_output_no_input_carries_names_the_route():
    with pytest.raises(ScientificToolchainContractError) as refused:
        _claim_node(("dg_sulfone_mecn",), ("dg_gauche_anti_sulfone_mecn",))
    assert refused.value.cause == "claim_output_unrendered"
    assert "rename the claim node's input_id" in (
        refused.value.next_legal_route
    )
