"""A completion receipt has an id space for what was found and not asked.

It had one for what was asked and missed (limitations) and none for the
other direction, so a host-detected anomaly could ride no receipt a
settlement reads. The new list is neither a finding nor a limitation:
a passed completion stays passed, the field is absent from the digest
when empty so every earlier receipt verifies, and each id carries the
signal and the standing status.
"""

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.analysis_completion import (
    AnalysisCompletionReceiptV1,
    _completion_receipt_body,
)


def _body(**extra):
    return {
        "schema_version": "chemsmart.analysis-completion-receipt.v1",
        "policy_sha256": "b" * 64,
        "task_spec_sha256": "c" * 64,
        "source_receipt_sha256s": ("d" * 64,),
        "status": "passed",
        "findings": (),
        **extra,
    }


@pytest.mark.capability("signal:stationary_point.unexpected_order")
def test_observations_ride_a_passed_completion_without_moving_its_digest():
    plain = AnalysisCompletionReceiptV1(
        **_body(), receipt_sha256=canonical_sha256(_body())
    )
    ids = ("anomaly:stationary_point.unexpected_order:unreplicated:0a1b2c3d",)
    with_observations = AnalysisCompletionReceiptV1(
        **_body(anomaly_output_ids=ids),
        receipt_sha256=canonical_sha256(_body(anomaly_output_ids=ids)),
    )
    assert with_observations.status == "passed"
    assert with_observations.limitation_output_ids == ()
    assert "anomaly_output_ids" not in _completion_receipt_body(plain)
    assert (
        _completion_receipt_body(with_observations)["anomaly_output_ids"]
        == ids
    )
    with pytest.raises(ContractError, match="anomaly output ids"):
        AnalysisCompletionReceiptV1(
            **_body(anomaly_output_ids=("b", "a")),
            receipt_sha256=canonical_sha256(
                _body(anomaly_output_ids=("b", "a"))
            ),
        )
