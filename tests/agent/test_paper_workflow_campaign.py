"""Focused paper-campaign measurement contracts."""

import pytest
from types import SimpleNamespace

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.analysis_claims import (
    AnalysisReportedQuantityV1,
    build_analysis_claim_record,
)
from chemsmart.agent.experiments.paper_workflow_campaign import (
    analysis_claim_record_from_live_result,
    PaperWorkflowObservationV1,
    WorkflowDependencyAnswerV1,
    WorkflowProjectAnswerV1,
    WorkflowSemanticNodeV1,
    build_paper_answer_key,
    build_paper_black_box_case,
    build_paper_workflow_answer_key,
    build_reference_quantity,
    grade_paper_numerical_claims,
    grade_paper_workflow_answer,
    PaperWorkflowOracleV1,
    grade_paper_workflow_result,
    prepare_benchmark_dispatch,
)
from chemsmart.agent.projects import project_document
from chemsmart.analysis.quantity_expressions import normalize_numeric_value
from chemsmart.analysis.result_quantities import ENERGY, make_quantity_value


def _case(*, case_id="paper-case", system_scale="small_molecule"):
    return build_paper_black_box_case(
        case_id=case_id,
        paper_doi="10.1000/example",
        article_title="Example calculation",
        source_locator="Methods, page 2",
        methods_excerpt="Compute a state-energy difference.",
        coordinate_sha256="a" * 64,
        coordinate_description="Two exact official XYZ artifacts.",
        requested_result="Report the spin gap in kcal/mol.",
        execution_policy="Run only supported approved nodes.",
        system_scale=system_scale,
    )


def _key(case, *, value=2.50, tolerance=0.10):
    reference = build_reference_quantity(
        quantity_id="spin-gap",
        expected_value=value,
        unit="kcal/mol",
        absolute_tolerance=tolerance,
        evidence_locator="private validated reference receipt",
    )
    return build_paper_answer_key(
        answer_key_id=f"{case.case_id}-answers",
        case=case,
        answer_source="host-prepared result withheld from the model",
        quantities=(reference,),
    )


def _claim_record(*, value, unit, claim_id="spin-gap", quantity_id="spin-gap"):
    canonical_value, canonical_unit, dimension = normalize_numeric_value(value, unit)
    quantity = make_quantity_value(
        quantity_id=quantity_id,
        source_value=value,
        source_unit=unit,
        value=canonical_value,
        unit=canonical_unit,
        dimension=dimension,
        evidence_ref="expression:spin-gap#" + "b" * 64,
    )
    claim = AnalysisReportedQuantityV1(
        claim_id=claim_id,
        source_kind="quantity_expression",
        source_receipt_sha256="c" * 64,
        quantity_id=quantity.quantity_id,
        quantity_value_sha256=quantity.value_sha256,
        display_value=value,
        display_unit=unit,
        canonical_value=quantity.value,
        canonical_unit=quantity.unit,
        dimension=ENERGY,
        data_kind="scalar",
    )
    return build_analysis_claim_record(
        task_spec_sha256="d" * 64, claims=(claim,)
    )


def test_scientific_term_grade_is_separator_insensitive():
    result = SimpleNamespace(
        public_transcript=(
            {
                "content": (
                    "Use def2SVP and extract oscillator strengths plus "
                    "excitation energies. RI CC2 is blocked_unsupported."
                ),
                "tool_calls": (
                    {"function": {"name": "plan_scientific_workflow"}},
                ),
            },
        ),
        final_text="PBE0/def2-SVP and RI-CC2 are preserved.",
        terminal_state="planned",
    )
    oracle = PaperWorkflowOracleV1(
        oracle_id="separator-normalization",
        required_tools=("plan_scientific_workflow",),
        required_transcript_terms=(
            "def2-svp",
            "excitation_energies",
            "oscillator_strengths",
        ),
        required_final_terms=("ri-cc2",),
        unsupported_terms=("ri-cc2",),
        require_preview=False,
    )
    grade = grade_paper_workflow_result(result=result, oracle=oracle)
    assert grade.strict_pass
    assert grade.evaluation_role == "diagnostic_structure_only"


def test_benchmark_dispatch_requires_a_matching_prepared_answer_key():
    case = _case()
    key = _key(case)
    task, eligibility = prepare_benchmark_dispatch(case=case, answer_key=key)

    assert eligibility.status == "eligible_answer_key_registered"
    assert eligibility.required_quantity_ids == ("spin-gap",)
    assert "2.5" not in task
    assert "private validated reference receipt" not in task
    assert case.article_title not in task
    assert case.paper_doi not in task
    assert case.source_locator not in task
    assert case.case_id not in task
    assert case.methods_excerpt in task
    assert case.coordinate_description in task
    assert case.requested_result in task
    assert case.execution_policy in task

    other_case = _case(case_id="other-case")
    with pytest.raises(ContractError, match="another black-box case"):
        prepare_benchmark_dispatch(case=other_case, answer_key=key)

    larger = _case(
        case_id="larger-numerical", system_scale="intermediate_or_larger"
    )
    larger_key = _key(larger)
    with pytest.raises(ContractError, match="restricted to small molecules"):
        prepare_benchmark_dispatch(case=larger, answer_key=larger_key)


def test_primary_grade_requires_the_numerical_claim_not_a_structural_pass():
    case = _case()
    key = _key(case)

    missing = grade_paper_numerical_claims(
        case=case, answer_key=key, claim_record=None
    )
    assert not missing.primary_pass
    assert missing.quantity_grades[0].status == "missing_claim"


def test_primary_grade_compares_compatible_units_with_declared_tolerance():
    case = _case()
    key = _key(case, value=2.50, tolerance=0.10)

    within = grade_paper_numerical_claims(
        case=case,
        answer_key=key,
        claim_record=_claim_record(value=10.50, unit="kJ/mol"),
    )
    outside = grade_paper_numerical_claims(
        case=case,
        answer_key=key,
        claim_record=_claim_record(value=12.0, unit="kJ/mol"),
    )

    assert within.primary_pass
    assert within.quantity_grades[0].status == "within_tolerance"
    assert not outside.primary_pass
    assert outside.quantity_grades[0].status == "outside_tolerance"


def test_primary_grade_uses_canonical_quantity_not_free_form_claim_label():
    case = _case()
    key = _key(case, value=2.50, tolerance=0.10)

    grade = grade_paper_numerical_claims(
        case=case,
        answer_key=key,
        claim_record=_claim_record(
            value=10.50,
            unit="kJ/mol",
            claim_id="reported-state-energy-difference",
            quantity_id="spin-gap",
        ),
    )

    assert grade.primary_pass


def test_primary_grade_matches_conservative_semantic_claim_tokens():
    case = _case()
    key = _key(case, value=2.50, tolerance=0.10)

    grade = grade_paper_numerical_claims(
        case=case,
        answer_key=key,
        claim_record=_claim_record(
            value=10.50,
            unit="kJ/mol",
            claim_id="reported-molar-spin-gap",
            quantity_id="gap-display",
        ),
    )

    assert grade.primary_pass


def test_primary_grade_uses_preregistered_identifier_without_paper_alias_code():
    case = _case()
    reference = build_reference_quantity(
        quantity_id="dimethyl-ether-energy",
        expected_value=-11.388320847866,
        unit="hartree",
        absolute_tolerance=1.0e-9,
        evidence_locator="private validated reference receipt",
        accepted_identifiers=("dme-sp-energy",),
    )
    key = build_paper_answer_key(
        answer_key_id="paper-case-alias-answers",
        case=case,
        answer_source="host-prepared result withheld from the model",
        quantities=(reference,),
    )

    grade = grade_paper_numerical_claims(
        case=case,
        answer_key=key,
        claim_record=_claim_record(
            value=-11.388320847866,
            unit="hartree",
            claim_id="dme_sp_energy_claim",
            quantity_id="dme_sp_energy",
        ),
    )

    assert grade.primary_pass


def test_answer_key_rejects_identifier_shared_by_two_quantities():
    case = _case()
    references = tuple(
        build_reference_quantity(
            quantity_id=quantity_id,
            expected_value=value,
            unit="hartree",
            absolute_tolerance=1.0e-9,
            evidence_locator="private validated reference receipt",
            accepted_identifiers=("shared-energy",),
        )
        for quantity_id, value in (("energy-a", -1.0), ("energy-b", -2.0))
    )

    with pytest.raises(ContractError, match="identify only one"):
        build_paper_answer_key(
            answer_key_id="paper-case-shared-alias",
            case=case,
            answer_source="host-prepared result withheld from the model",
            quantities=references,
        )


def test_primary_grade_collapses_equivalent_multi_unit_claim_presentations():
    case = _case()
    key = _key(case, value=2.50, tolerance=0.10)
    canonical_value, canonical_unit, dimension = normalize_numeric_value(
        2.50, "kcal/mol"
    )
    quantity = make_quantity_value(
        quantity_id="delta-energy",
        source_value=2.50,
        source_unit="kcal/mol",
        value=canonical_value,
        unit=canonical_unit,
        dimension=dimension,
        evidence_ref="expression:spin-gap#" + "b" * 64,
    )
    claims = tuple(
        AnalysisReportedQuantityV1(
            claim_id=claim_id,
            source_kind="quantity_expression",
            source_receipt_sha256="c" * 64,
            quantity_id=quantity.quantity_id,
            quantity_value_sha256=quantity.value_sha256,
            display_value=value,
            display_unit=unit,
            canonical_value=quantity.value,
            canonical_unit=quantity.unit,
            dimension=ENERGY,
            data_kind="scalar",
        )
        for claim_id, value, unit in (
            ("spin-gap-kcal", 2.50, "kcal/mol"),
            ("spin-gap-kj", 10.46, "kJ/mol"),
        )
    )
    record = build_analysis_claim_record(
        task_spec_sha256="d" * 64, claims=claims
    )

    grade = grade_paper_numerical_claims(
        case=case, answer_key=key, claim_record=record
    )

    assert grade.primary_pass


def test_primary_grade_normalizes_hartree_identifier_abbreviation():
    case = _case()
    reference = build_reference_quantity(
        quantity_id="isomerization-energy",
        expected_value=0.003,
        unit="hartree",
        absolute_tolerance=1.0e-9,
        evidence_locator="private validated reference receipt",
        accepted_identifiers=("delta-e-hartree",),
    )
    key = build_paper_answer_key(
        answer_key_id="paper-case-hartree-alias",
        case=case,
        answer_source="host-prepared result withheld from the model",
        quantities=(reference,),
    )

    grade = grade_paper_numerical_claims(
        case=case,
        answer_key=key,
        claim_record=_claim_record(
            value=0.003,
            unit="hartree",
            claim_id="isomerization-delta-e-eh",
            quantity_id="delta-e",
        ),
    )

    assert grade.primary_pass


def test_live_primary_grade_uses_host_claim_record_not_model_prose():
    record = _claim_record(value=10.50, unit="kJ/mol")
    envelope = {
        "schema_version": "chemsmart.tool-result.v1",
        "status": "ok",
        "tool": "record_analysis_claims",
        "result": {
            "schema_version": record.schema_version,
            "task_spec_sha256": record.task_spec_sha256,
            "claims": [claim.__dict__ for claim in record.claims],
            "status": record.status,
            "receipt_sha256": record.receipt_sha256,
        },
    }
    result = SimpleNamespace(
        task_spec_sha256=record.task_spec_sha256,
        public_transcript=(
            {"role": "assistant", "content": "The answer is 999 kcal/mol."},
            {"role": "tool", "content": __import__("json").dumps(envelope)},
        ),
    )

    observed = analysis_claim_record_from_live_result(result)

    assert observed == record


def _orca_workflow(*, opt_node_id, analysis_node_id, jobtype="opt"):
    project = WorkflowProjectAnswerV1(
        project_id=f"{opt_node_id}-project",
        document=project_document(
            program="orca",
            sections={
                "gas": {
                    "method": "B3LYP",
                    "basis": "def2-SVP",
                    "freq": True,
                }
            },
        ),
    )
    opt = WorkflowSemanticNodeV1(
        node_id=opt_node_id,
        node_kind="calculation",
        program="orca",
        jobtype=jobtype,
        project_id=project.project_id,
        input_geometry_sha256s=("a" * 64,),
        output_semantics=("optimized_geometry", "vibrational_frequencies"),
    )
    analysis = WorkflowSemanticNodeV1(
        node_id=analysis_node_id,
        node_kind="analysis",
        operation="result_extraction",
        dependencies=(
            WorkflowDependencyAnswerV1(
                source_node_id=opt_node_id,
                source_output="vibrational_frequencies",
                target_input="frequencies",
            ),
        ),
        output_semantics=("frequency_table",),
    )
    return (project,), tuple(sorted((opt, analysis), key=lambda item: item.node_id))


def test_intermediate_case_uses_exact_program_relative_yaml_cli_dag_answer():
    case = _case(
        case_id="orca-intermediate", system_scale="intermediate_or_larger"
    )
    projects, nodes = _orca_workflow(
        opt_node_id="reference-opt-freq", analysis_node_id="reference-extract"
    )
    key = build_paper_workflow_answer_key(
        answer_key_id="orca-intermediate-answer",
        case=case,
        answer_source="expert-prepared ChemSmart ORCA project and DAG",
        projects=projects,
        nodes=nodes,
    )
    task, eligibility = prepare_benchmark_dispatch(case=case, answer_key=key)
    assert task
    assert eligibility.evaluation_mode == "canonical_yaml_cli_dag"

    renamed_projects, renamed_nodes = _orca_workflow(
        opt_node_id="model-optimization", analysis_node_id="model-parser"
    )
    correct = grade_paper_workflow_answer(
        case=case,
        answer_key=key,
        observation=PaperWorkflowObservationV1(
            schema_version="chemsmart.paper-workflow-observation.v1",
            projects=renamed_projects,
            nodes=renamed_nodes,
        ),
    )
    wrong_projects, wrong_nodes = _orca_workflow(
        opt_node_id="model-hessian", analysis_node_id="model-parser", jobtype="hess"
    )
    wrong = grade_paper_workflow_answer(
        case=case,
        answer_key=key,
        observation=PaperWorkflowObservationV1(
            schema_version="chemsmart.paper-workflow-observation.v1",
            projects=wrong_projects,
            nodes=wrong_nodes,
        ),
    )

    assert correct.strict_pass
    assert not wrong.strict_pass
    assert wrong.missing_node_signatures
    assert wrong.unexpected_node_signatures


def test_intermediate_strict_grade_distinguishes_cbs_law_parameters():
    case = _case(
        case_id="cbs-intermediate", system_scale="intermediate_or_larger"
    )
    project = WorkflowProjectAnswerV1(
        project_id="orca-qz",
        document=project_document(
            program="orca",
            sections={"solv": {"ab_initio": "dlpno-ccsd(t1)", "basis": "cc-pvqz-dk"}},
        ),
    )
    calculation = WorkflowSemanticNodeV1(
        node_id="sp-qz",
        node_kind="calculation",
        program="orca",
        jobtype="sp",
        project_id=project.project_id,
        input_geometry_sha256s=("a" * 64,),
        output_semantics=("program_result",),
    )
    reference_analysis = WorkflowSemanticNodeV1(
        node_id="scf-cbs",
        node_kind="analysis",
        operation="quantity_expression",
        dependencies=(
            WorkflowDependencyAnswerV1("sp-qz", "program_result", "program_result"),
        ),
        output_semantics=("scf_energy",),
        semantic_parameters=(
            ("cardinal_numbers", (4, 5)),
            ("extrapolation_exponent", 3.9),
            ("operation", "scf_inverse_power_cbs_limit"),
        ),
    )
    key = build_paper_workflow_answer_key(
        answer_key_id="cbs-intermediate-answer",
        case=case,
        answer_source="expert-prepared inverse-power CBS contract",
        projects=(project,),
        nodes=(calculation, reference_analysis),
    )
    wrong_analysis = WorkflowSemanticNodeV1(
        node_id="renamed-cbs",
        node_kind="analysis",
        operation="quantity_expression",
        dependencies=(
            WorkflowDependencyAnswerV1("renamed-sp", "program_result", "program_result"),
        ),
        output_semantics=("scf_energy",),
        semantic_parameters=(
            ("cardinal_numbers", (3, 4)),
            ("extrapolation_exponent", 3.0),
            ("operation", "correlation_inverse_power_cbs_limit"),
        ),
    )
    renamed_calculation = WorkflowSemanticNodeV1(
        node_id="renamed-sp",
        node_kind="calculation",
        program="orca",
        jobtype="sp",
        project_id=project.project_id,
        input_geometry_sha256s=("a" * 64,),
        output_semantics=("program_result",),
    )

    grade = grade_paper_workflow_answer(
        case=case,
        answer_key=key,
        observation=PaperWorkflowObservationV1(
            schema_version="chemsmart.paper-workflow-observation.v1",
            projects=(project,),
            nodes=(renamed_calculation, wrong_analysis),
        ),
    )

    assert not grade.strict_pass
    assert grade.missing_node_signatures
    assert grade.unexpected_node_signatures


def test_orca_sp_grading_compares_effective_settings_not_section_spelling():
    case = _case(
        case_id="orca-sp-sections", system_scale="intermediate_or_larger"
    )
    legacy = WorkflowProjectAnswerV1(
        project_id="reference-sp",
        document=project_document(
            program="orca",
            sections={
                "solv": {
                    "ab_initio": "dlpno-ccsd(t1)",
                    "basis": "cc-pvtz-dk",
                }
            },
        ),
    )
    reference_node = WorkflowSemanticNodeV1(
        node_id="reference-node",
        node_kind="calculation",
        program="orca",
        jobtype="sp",
        project_id=legacy.project_id,
        input_geometry_sha256s=("a" * 64,),
        output_semantics=("program_result",),
    )
    key = build_paper_workflow_answer_key(
        answer_key_id="orca-sp-section-answer",
        case=case,
        answer_source="expert-prepared effective ORCA SP settings",
        projects=(legacy,),
        nodes=(reference_node,),
    )

    explicit_stage = WorkflowProjectAnswerV1(
        project_id="model-sp",
        document=project_document(
            program="orca",
            sections={
                "gas": {},
                "sp": {
                    "ab_initio": "DLPNO-CCSD(T1)",
                    "basis": "cc-pVTZ-DK",
                },
            },
        ),
    )
    model_node = WorkflowSemanticNodeV1(
        node_id="model-node",
        node_kind="calculation",
        program="orca",
        jobtype="sp",
        project_id=explicit_stage.project_id,
        input_geometry_sha256s=("a" * 64,),
        output_semantics=("program_result",),
    )
    equivalent = grade_paper_workflow_answer(
        case=case,
        answer_key=key,
        observation=PaperWorkflowObservationV1(
            schema_version="chemsmart.paper-workflow-observation.v1",
            projects=(explicit_stage,),
            nodes=(model_node,),
        ),
    )
    phase_keyed = WorkflowProjectAnswerV1(
        project_id="phase-keyed-sp",
        document=project_document(
            program="orca",
            sections={
                "gas": {
                    "ab_initio": "dlpno-ccsd(t1)",
                    "basis": "cc-pvtz-dk",
                }
            },
        ),
    )
    phase_keyed_node = WorkflowSemanticNodeV1(
        node_id="phase-keyed-node",
        node_kind="calculation",
        program="orca",
        jobtype="sp",
        project_id=phase_keyed.project_id,
        input_geometry_sha256s=("a" * 64,),
        output_semantics=("program_result",),
    )
    phase_equivalent = grade_paper_workflow_answer(
        case=case,
        answer_key=key,
        observation=PaperWorkflowObservationV1(
            schema_version="chemsmart.paper-workflow-observation.v1",
            projects=(phase_keyed,),
            nodes=(phase_keyed_node,),
        ),
    )

    assert equivalent.strict_pass
    assert phase_equivalent.strict_pass
