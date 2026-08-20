"""A waiting node must say what it waits for, and a producer what it feeds.

Observed live: every waiting node in a real DeepSeek workflow received the
identical string ``await producer outputs``.  Six nodes, six copies, no
producer named and no output named -- while the host held the exact edges the
model had just submitted.  A model that has to recall its own DAG re-plans,
and re-planning is where the stage carrying the answer gets dropped.

These tests pin that the frontier states the branching it already knows.
"""

import pytest

from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    build_scientific_toolchain_plan,
    project_scientific_toolchain_frontier,
)
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)


def _extraction(node_id, producer, output_id):
    return AnalysisNodeIntentV1(
        node_id=node_id,
        analysis_kind="result_extraction",
        dependencies=(producer,),
        inputs=(
            AnalysisInputIntentV1(
                input_id=f"raw-{output_id}",
                source_kind="program_output",
                producer_node_id=producer,
                producer_output_id=f"result-{output_id}",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id=output_id, selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id=output_id, quantity_kind="energy", unit="hartree"
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )


def _limit(node_id, sources):
    return AnalysisNodeIntentV1(
        node_id=node_id,
        analysis_kind="quantity_expression",
        dependencies=tuple(sorted(producer for producer, _ in sources)),
        inputs=tuple(
            AnalysisInputIntentV1(
                input_id=output_id,
                source_kind="analysis_output",
                producer_node_id=producer,
                producer_output_id=output_id,
            )
            for producer, output_id in sorted(sources, key=lambda s: s[1])
        ),
        selectors=(),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e_limit", quantity_kind="energy", unit="hartree"
            ),
        ),
        expression_nodes=(
            {
                "node_id": "e_limit",
                "operation": "exponential_cbs_limit",
                "input_ids": sorted(output_id for _, output_id in sources),
            },
        ),
        expression_output_node_ids=("e_limit",),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )


def _calculation(node_id, output_id):
    return CommandNodeIntentV1(
        node_id=node_id,
        program="pyscf",
        jobtype="sp",
        project_role="sp",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geometry",
                artifact_class="geometry_xyz",
                artifact_id="start.xyz",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id=output_id, artifact_class="pyscf_hdf5"
            ),
        ),
        unresolved_fields=(),
    )


def _plan(analysis_nodes, calculations, observables):
    return build_scientific_toolchain_plan(
        plan_id="probe-plan",
        workflow_id="probe",
        command_workflow_draft_sha256="a" * 64,
        calculation_nodes=tuple(
            _calculation(node_id, outputs[0])
            for node_id, outputs in observables
        ),
        calculation_observables=dict(observables),
        analysis_nodes=analysis_nodes,
        required_output_ids=("e_limit",),
    )


@pytest.fixture
def frontier():
    calculations = ("sp-dz", "sp-qz", "sp-tz")
    observables = tuple(
        (node, (f"result-e_{node[-2:]}",)) for node in calculations
    )
    analysis = (
        _extraction("extract-dz", "sp-dz", "e_dz"),
        _extraction("extract-tz", "sp-tz", "e_tz"),
        _extraction("extract-qz", "sp-qz", "e_qz"),
        _limit(
            "cbs",
            (
                ("extract-dz", "e_dz"),
                ("extract-tz", "e_tz"),
                ("extract-qz", "e_qz"),
            ),
        ),
    )
    return project_scientific_toolchain_frontier(
        _plan(analysis, calculations, observables),
        actionable_calculation_node_ids=calculations,
    )


def _node(frontier, node_id):
    return next(
        item for item in frontier["nodes"] if item["node_id"] == node_id
    )


def test_a_waiting_node_names_the_producer_output_it_needs(frontier):
    node = _node(frontier, "extract-dz")
    assert node["state"] == "waiting_for_artifact"
    assert node["waiting_on"] == ("sp-dz",)
    assert node["unsatisfied_inputs"][0]["producer_output_id"] == "result-e_dz"
    assert "sp-dz.result-e_dz" in node["next_action"]


def test_a_fan_in_node_names_every_producer_not_just_one(frontier):
    node = _node(frontier, "cbs")
    assert node["waiting_on"] == ("extract-dz", "extract-qz", "extract-tz")
    named = node["next_action"]
    for producer in ("extract-dz.e_dz", "extract-qz.e_qz", "extract-tz.e_tz"):
        assert producer in named, named


def test_branching_is_visible_from_the_producer(frontier):
    """Fan-out is a DAG property the model should not have to recall."""

    assert _node(frontier, "sp-dz")["feeds"] == ("extract-dz",)
    assert _node(frontier, "extract-dz")["feeds"] == ("cbs",)
    assert _node(frontier, "cbs")["feeds"] == ()


def test_no_two_waiting_nodes_share_one_contentless_sentence(frontier):
    """The defect this replaces: one identical string for every waiting node."""

    reasons = [
        item["next_action"]
        for item in frontier["nodes"]
        if item["state"] == "waiting_for_artifact"
    ]
    assert len(reasons) == len(set(reasons)), reasons
    assert "await producer outputs" not in reasons


def test_completed_calculation_is_not_offered_for_rerun():
    calculations = ("sp-dz",)
    observables = (("sp-dz", ("result-e_dz",)),)
    analysis = (
        _extraction("extract-dz", "sp-dz", "e_dz"),
        _limit("cbs", (("extract-dz", "e_dz"),)),
    )
    projected = project_scientific_toolchain_frontier(
        _plan(analysis, calculations, observables),
        completed_calculation_node_ids=calculations,
    )

    assert _node(projected, "sp-dz")["state"] == "completed"
    assert projected["completed_node_ids"] == ("sp-dz",)
    # The completed calculation is not re-offered; what becomes actionable is
    # the extraction its completion just satisfied.
    assert "sp-dz" not in projected["actionable_node_ids"]
    assert projected["actionable_node_ids"] == ("extract-dz",)


def test_a_blocked_upstream_node_names_which_parent_is_blocked():
    calculations = ("sp-dz",)
    observables = (("sp-dz", ("result-e_dz",)),)
    unsupported = AnalysisNodeIntentV1(
        node_id="unsupported",
        analysis_kind="scientific_validation",
        dependencies=("sp-dz",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp-dz",
                producer_output_id="result-e_dz",
            ),
        ),
        selectors=(),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="verdict", quantity_kind="flag", unit="1"
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="blocked_unsupported",
        blocked_reason="no local-mode analysis exists",
    )
    downstream = AnalysisNodeIntentV1(
        node_id="report",
        analysis_kind="claim_rendering",
        dependencies=("unsupported",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="v",
                source_kind="analysis_output",
                producer_node_id="unsupported",
                producer_output_id="verdict",
            ),
        ),
        selectors=(),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e_limit", quantity_kind="record", unit="1"
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )
    frontier = project_scientific_toolchain_frontier(
        _plan((unsupported, downstream), calculations, observables),
        actionable_calculation_node_ids=calculations,
    )
    node = _node(frontier, "report")
    assert node["state"] == "blocked_upstream"
    assert node["waiting_on"] == ("unsupported",)
    assert "unsupported is blocked" in node["next_action"]


def test_the_frontier_is_derived_never_asserted_by_the_model(frontier):
    """Every edge reported here came from the submitted plan."""

    submitted = {
        "sp-dz",
        "sp-qz",
        "sp-tz",
        "extract-dz",
        "extract-qz",
        "extract-tz",
        "cbs",
    }
    for item in frontier["nodes"]:
        assert set(item.get("waiting_on", ())) <= submitted
        assert set(item["feeds"]) <= submitted
