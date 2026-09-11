"""What the host could state about a node reaches the human who decides.

NOVEL-3 ino3 (2026-09-05): ORCA's 3N geometry-iteration cap was stated
62 times in one goal's compile results and reached no review, no bundle
and no human; two nodes then stopped at the cap. In the same goal a
quartet project promoted with ``geom_maxiter: 150`` was re-promoted as
a ``-v2`` without it after an unrelated refusal, and nothing said the
lever had been lost.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.execution import (
    approve_workflow_execution_review,
    build_workflow_execution_review,
)
from chemsmart.agent.tui.review import render_review_blocks
from tests.agent.test_exact_execution_approval_chain import _review

_THREE_N = (
    "ORCA's default geometry-iteration cap for this 21-atom input is "
    "3N = 63; a floppy or metal-centred system may need more, and "
    "geom_maxiter on the project raises it"
)


def _rebuilt(review, node_observations):
    return build_workflow_execution_review(
        request=review.request,
        scientific_plan=review.scientific_plan,
        materialized_workflow=review.materialized_workflow,
        execution_resources=review.execution_resources,
        execution_envelope=review.execution_envelope,
        environment_bindings=review.environment_bindings,
        node_reviews=review.node_reviews,
        stationary_point_policy=review.stationary_point_policy,
        node_observations=node_observations,
    )


@pytest.mark.capability("rule:review.displays_host_observations")
def test_node_observations_ride_the_review_the_display_and_the_bundle(
    tmp_path,
):
    plain = _review(tmp_path)
    assert plain.node_observations == ()
    # An empty observation list keeps the original canonical body.
    assert (
        _rebuilt(plain, [{"node_id": "sp", "observations": ()}]).review_sha256
        == plain.review_sha256
    )

    stated = _rebuilt(plain, [{"node_id": "sp", "observations": (_THREE_N,)}])
    assert stated.review_sha256 != plain.review_sha256
    assert stated.node_observations[0]["observations"] == [_THREE_N]

    panels = [
        block
        for block in render_review_blocks(stated)
        if getattr(block, "title", None) == "sp · host observations"
    ]
    assert len(panels) == 1
    assert "3N = 63" in str(panels[0].renderable)

    bundle = approve_workflow_execution_review(
        stated,
        approval_id="approval-sp",
        approved_review_sha256=stated.review_sha256,
        actor="human",
        resolution_id="resolution-sp",
    )
    assert bundle.node_observations == stated.node_observations


def test_a_re_promotion_that_drops_a_field_is_told_so(tmp_path):
    """ino3's shape: v1 carried geom_maxiter, v2 did not."""

    from tests.agent.test_a_guide_opens_when_something_asks import _host

    host = _host(tmp_path)
    first = host._render_project_yaml(
        "t1",
        {
            "program": "orca",
            "sections": {
                "gas": {"functional": "b3lyp", "basis": "def2-svp"},
                "opt": {"geom_maxiter": 150},
            },
        },
    )
    promoted = host._promote_project_yaml(
        "t1",
        {
            "render_receipt_sha256": first.receipt_sha256,
            "artifact_id": "quartet-project",
        },
    )
    assert promoted["observations"] == ()

    second = host._render_project_yaml(
        "t2",
        {
            "program": "orca",
            "sections": {
                "gas": {"functional": "b3lyp", "basis": "def2-svp"},
                "opt": {"opt_convergence": "tight"},
            },
        },
    )
    again = host._promote_project_yaml(
        "t2",
        {
            "render_receipt_sha256": second.receipt_sha256,
            "artifact_id": "quartet-project-v2",
        },
    )
    (observation,) = again["observations"]
    assert observation.startswith("against 'quartet-project' (orca): ")
    assert "opt.geom_maxiter (150)" in observation
    assert "restates it" in observation

    # A project of another program, or a section only one carries, says
    # nothing: the sentence is for the same role.
    third = host._render_project_yaml(
        "t3",
        {
            "program": "orca",
            "sections": {"gas": {"functional": "b3lyp", "basis": "def2-svp"}},
        },
    )
    sp_only = host._promote_project_yaml(
        "t3",
        {
            "render_receipt_sha256": third.receipt_sha256,
            "artifact_id": "sp-project",
        },
    )
    assert sp_only["observations"] == ()
