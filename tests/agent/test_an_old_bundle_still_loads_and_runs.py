"""Every packet approved before the analysis chain existed keeps its bytes.

The toolchain plan is an additive review/bundle field. A calculation-only
review built today must produce byte-identical canonical bodies to one
built before the field existed, and packets serialized without the key
must load and re-verify their original digests.
"""

from __future__ import annotations

import json

from chemsmart.agent._contracts import canonical_data
from chemsmart.agent.execution import (
    approve_workflow_execution_review,
    workflow_execution_review_json,
)
from chemsmart.agent.live_session import (
    load_workflow_execution_review,
)
from tests.agent.test_exact_execution_approval_chain import _review


def test_a_chainless_packet_has_no_chain_key_and_round_trips(tmp_path):
    review = _review(tmp_path)

    body = canonical_data(review)
    assert "scientific_toolchain_plan" not in json.dumps(body) or (
        body.get("scientific_toolchain_plan") is None
    )
    assert "extract" not in workflow_execution_review_json(review)

    path = tmp_path / "review.json"
    path.write_text(workflow_execution_review_json(review), encoding="utf-8")
    assert load_workflow_execution_review(path) == review

    bundle = approve_workflow_execution_review(
        review,
        approval_id="approval-water-sp",
        approved_review_sha256=review.review_sha256,
        actor="human",
        resolution_id="resolution-water-sp",
    )
    assert bundle.scientific_toolchain_plan is None
    assert (
        bundle.frozen_workflow_approval.scientific_toolchain_plan_sha256 == ""
    )
