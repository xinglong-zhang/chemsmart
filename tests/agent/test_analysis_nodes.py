"""Focused contracts for the program-neutral scientific analysis overlay."""

from dataclasses import replace

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.analysis_nodes import (
    AnalysisContractError,
    AnalysisInputRefV1,
    AnalysisOutputSpecV1,
    PySCFResultParserAdapterV1,
    ResultParserRegistryV1,
    ResultQuantitySelectorV1,
    build_analysis_node_spec,
    build_default_result_parser_registry,
    build_scientific_analysis_plan,
    execute_result_extraction_node,
)
from chemsmart.analysis.result_quantities import (
    DIMENSIONLESS,
    ENERGY,
    QuantityExtractionReceiptV1,
    QuantitySelectorV1,
    canonical_quantity_sha256,
    make_quantity_value,
)


_TASK = "a" * 64
_WORKFLOW_PLAN = "b" * 64
_SOURCE = "c" * 64
_SOURCE_RECEIPT = "d" * 64
_ANALYSIS_PLAN = "e" * 64


def _selector(quantity_id="energy", selector="energy"):
    return ResultQuantitySelectorV1(
        schema_version="chemsmart.result-quantity-selector.v1",
        quantity_id=quantity_id,
        selector=selector,
    )


def _program_input(*, source_sha256=_SOURCE):
    return AnalysisInputRefV1(
        input_id="result",
        source_kind="program_artifact",
        producer_node_id="opt",
        producer_output_id="structured_result",
        source_id="result.opt.hdf5",
        source_sha256=source_sha256,
        source_receipt_sha256=_SOURCE_RECEIPT,
        program="pyscf",
    )


def _output(quantity_id="energy", unit="hartree", dimension=ENERGY):
    return AnalysisOutputSpecV1(
        quantity_id=quantity_id,
        unit=unit,
        dimension=dimension,
    )


def _extraction_node(
    *,
    selector=None,
    output=None,
    source_sha256=_SOURCE,
    support_state="resolvable",
    blocked_reason="",
):
    return build_analysis_node_spec(
        node_id="extract",
        task_spec_sha256=_TASK,
        workflow_id="wf.probe",
        scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
        analysis_kind="result_extraction",
        inputs=(_program_input(source_sha256=source_sha256),),
        selectors=(selector or _selector(),),
        outputs=(output or _output(),),
        evidence_requirements=("validated-result",),
        support_state=support_state,
        blocked_reason=blocked_reason,
    )


def _analysis_input(producer, output_id, *, input_id="upstream"):
    return AnalysisInputRefV1(
        input_id=input_id,
        source_kind="analysis_quantity",
        producer_node_id=producer,
        producer_output_id=output_id,
        source_id=output_id,
        source_sha256="f" * 64,
        source_receipt_sha256="1" * 64,
    )


def _expression_node(node_id="derive", producer="extract", output_id="energy"):
    return build_analysis_node_spec(
        node_id=node_id,
        task_spec_sha256=_TASK,
        workflow_id="wf.probe",
        scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
        analysis_kind="quantity_expression",
        inputs=(_analysis_input(producer, output_id),),
        outputs=(_output("reported-energy"),),
        evidence_requirements=("dimensional-expression",),
    )


def _artifact(tmp_path, *, sha256=_SOURCE):
    path = tmp_path / "result.h5"
    path.write_bytes(b"not parsed by the patched adapter")
    return TrustedArtifactRefV1(
        artifact_id="result.opt.hdf5",
        kind="pyscf_hdf5",
        sha256=sha256,
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value="artifacts/result.h5",
    )


def _extraction_receipt(artifact):
    quantity = make_quantity_value(
        quantity_id="energy",
        source_value=-56.7,
        source_unit="hartree",
        value=-56.7,
        unit="hartree",
        dimension=ENERGY,
        evidence_ref=f"artifact:{artifact.artifact_id}#{artifact.sha256}",
    )
    body = {
        "schema_version": "chemsmart.quantity-extraction-receipt.v1",
        "artifact_id": artifact.artifact_id,
        "artifact_sha256": artifact.sha256,
        "program": "pyscf",
        "parser_id": "chemsmart.io.pyscf.output.PySCFOutput",
        "quantities": (quantity,),
        "status": "extracted",
    }
    return QuantityExtractionReceiptV1(
        **body, receipt_sha256=canonical_quantity_sha256(body)
    )


def test_selector_syntax_is_program_neutral_but_not_a_query_language():
    assert _selector("electronic-entropy", "electronic_entropy").selector == (
        "electronic_entropy"
    )
    with pytest.raises(AnalysisContractError, match="canonical public symbol"):
        _selector(selector="results/energies[-1]")
    with pytest.raises(AnalysisContractError, match="canonical public symbol"):
        _selector(selector="Energy")


def test_analysis_plan_is_bound_and_canonically_topological():
    extraction = _extraction_node()
    expression = _expression_node()
    plan = build_scientific_analysis_plan(
        analysis_plan_id="analysis.probe",
        task_spec_sha256=_TASK,
        workflow_id="wf.probe",
        scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
        nodes=(expression, extraction),
        required_output_quantity_ids=("reported-energy",),
    )
    assert tuple(node.node_id for node in plan.nodes) == ("extract", "derive")
    assert plan.required_output_quantity_ids == ("reported-energy",)
    assert len(plan.plan_sha256) == 64
    with pytest.raises(AnalysisContractError, match="digest mismatch"):
        replace(plan, plan_sha256="0" * 64)


def test_required_output_survives_as_a_blocked_analysis_node():
    node = _extraction_node(
        support_state="blocked_unsupported",
        blocked_reason="the selected parser does not expose this quantity",
    )
    plan = build_scientific_analysis_plan(
        analysis_plan_id="analysis.blocked",
        task_spec_sha256=_TASK,
        workflow_id="wf.probe",
        scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
        nodes=(node,),
        required_output_quantity_ids=("energy",),
    )
    assert plan.nodes[0].support_state == "blocked_unsupported"


def test_analysis_plan_rejects_cycles():
    left = build_analysis_node_spec(
        node_id="left",
        task_spec_sha256=_TASK,
        workflow_id="wf.probe",
        scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
        analysis_kind="quantity_expression",
        inputs=(_analysis_input("right", "right-value"),),
        outputs=(_output("left-value"),),
    )
    right = build_analysis_node_spec(
        node_id="right",
        task_spec_sha256=_TASK,
        workflow_id="wf.probe",
        scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
        analysis_kind="quantity_expression",
        inputs=(_analysis_input("left", "left-value"),),
        outputs=(_output("right-value"),),
    )
    with pytest.raises(AnalysisContractError, match="contains a cycle"):
        build_scientific_analysis_plan(
            analysis_plan_id="analysis.cycle",
            task_spec_sha256=_TASK,
            workflow_id="wf.probe",
            scientific_workflow_plan_sha256=_WORKFLOW_PLAN,
            nodes=(left, right),
        )


def test_pyscf_adapter_wraps_the_existing_provenance_reader(
    monkeypatch, tmp_path
):
    artifact = _artifact(tmp_path)
    expected = _extraction_receipt(artifact)
    observed = {}

    def fake_extract(*, request, artifact_path):
        observed["request"] = request
        observed["path"] = artifact_path
        return expected

    monkeypatch.setattr(
        "chemsmart.agent.analysis_nodes.extract_pyscf_quantities",
        fake_extract,
    )
    adapter = PySCFResultParserAdapterV1()
    receipt = adapter.extract(artifact=artifact, selectors=(_selector(),))
    assert receipt == expected
    assert observed["path"] == artifact.path
    assert isinstance(observed["request"].selectors[0], QuantitySelectorV1)


def test_registry_identity_is_order_independent():
    first = ResultParserRegistryV1((PySCFResultParserAdapterV1(),))
    second = ResultParserRegistryV1()
    second.register(PySCFResultParserAdapterV1())
    assert first.registry_sha256 == second.registry_sha256


def test_registered_extraction_produces_a_digest_bound_node_receipt(
    monkeypatch, tmp_path
):
    artifact = _artifact(tmp_path)
    expected = _extraction_receipt(artifact)
    monkeypatch.setattr(
        "chemsmart.agent.analysis_nodes.extract_pyscf_quantities",
        lambda **_: expected,
    )
    execution, extraction = execute_result_extraction_node(
        analysis_plan_sha256=_ANALYSIS_PLAN,
        node=_extraction_node(),
        artifact=artifact,
        registry=build_default_result_parser_registry(),
    )
    assert extraction == expected
    assert execution.status == "derived"
    assert execution.component_receipt_sha256s == (
        expected.receipt_sha256,
    )
    assert execution.outputs[0].value_sha256 == (
        expected.quantities[0].value_sha256
    )
    assert len(execution.receipt_sha256) == 64


def test_unregistered_selector_is_honestly_blocked(tmp_path):
    artifact = _artifact(tmp_path)
    selector = _selector("population", "population")
    node = _extraction_node(
        selector=selector,
        output=_output("population", "1", DIMENSIONLESS),
    )
    execution, extraction = execute_result_extraction_node(
        analysis_plan_sha256=_ANALYSIS_PLAN,
        node=node,
        artifact=artifact,
        registry=build_default_result_parser_registry(),
    )
    assert extraction is None
    assert execution.status == "blocked_unsupported"
    assert not execution.outputs
    assert execution.findings[0].rule_id == (
        "result-parser.selector-unsupported"
    )


def test_substituted_artifact_fails_before_parser_dispatch(tmp_path):
    artifact = _artifact(tmp_path, sha256="9" * 64)
    execution, extraction = execute_result_extraction_node(
        analysis_plan_sha256=_ANALYSIS_PLAN,
        node=_extraction_node(),
        artifact=artifact,
        registry=build_default_result_parser_registry(),
    )
    assert extraction is None
    assert execution.status == "failed"
    assert execution.findings[0].rule_id == (
        "result-parser.artifact-binding-mismatch"
    )


def test_parser_output_must_match_declared_units(monkeypatch, tmp_path):
    artifact = _artifact(tmp_path)
    expected = _extraction_receipt(artifact)
    monkeypatch.setattr(
        "chemsmart.agent.analysis_nodes.extract_pyscf_quantities",
        lambda **_: expected,
    )
    node = _extraction_node(output=_output("energy", "electronvolt", ENERGY))
    execution, extraction = execute_result_extraction_node(
        analysis_plan_sha256=_ANALYSIS_PLAN,
        node=node,
        artifact=artifact,
        registry=build_default_result_parser_registry(),
    )
    assert extraction == expected
    assert execution.status == "failed"
    assert execution.findings[0].rule_id == (
        "result-parser.output-contract-mismatch"
    )

