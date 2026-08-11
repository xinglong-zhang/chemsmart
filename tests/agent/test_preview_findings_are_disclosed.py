"""A refused preview must say what it objected to, not only that it did.

Measured on a live DeepSeek session (W23, 56 turns, 5.38M input tokens): the
model compiled Gaussian four times and ORCA twice, and every preview came back

    program_validation_status: "invalid"
    critical_finding_sha256s: ["94a63093...", "e3af4053...", ...]

with no text anywhere in the response and no tool able to resolve a digest. It
got PySCF right on its first attempt, gave up on the other two, and left the
dead branches in its plan -- which made the whole plan unfreezable, because
freezing requires every materialized node to be previewed.

The bodies were never missing. ``execute_safe_preview`` computed them in order
to hash them and dropped them on the next line. ``_public_validator_findings``
had been reading a ``findings`` attribute that no receipt ever set, so the
disclosure that existed rendered an empty list.

These gates pin the whole chain, because a break anywhere in it restores
silence.
"""

import inspect

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.preflight import (
    PreviewRuleFindingV1,
    validator_receipt_from_safe_preview,
)
from chemsmart.agent.preview import SafePreviewReceiptV1
from chemsmart.agent.program_verifiers import PreviewValidationFindingV1
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _public_validator_findings,
)

_FINDING = PreviewValidationFindingV1(
    rule_id="gaussian.route.functional_mismatch",
    field="route_functional",
    expected="hf",
    observed="b3lyp",
    evidence_ref="water_opt_fake.com",
)


def _failed_preview(**overrides) -> SafePreviewReceiptV1:
    body = {
        "schema_version": "chemsmart.safe-preview-receipt.v1",
        "invocation_sha256": "1" * 64,
        "observed_argv_sha256": "2" * 64,
        "input_sha256": "3" * 64,
        "project_sha256": "4" * 64,
        "exit_status": 0,
        "fake_mode": True,
        "no_scratch_mode": True,
        "scheduler_test_mode": False,
        "artifacts": (),
        "artifact_set_sha256": "5" * 64,
        "output_sha256": "6" * 64,
        "exception_class": "",
        "program_validation_receipt_sha256": "7" * 64,
        "program_validation_status": "invalid",
        "critical_finding_sha256s": (canonical_sha256(_FINDING),),
        "status": "failed",
        "rule_ids": ("preview.program_validator.red",),
    }
    body.update(overrides)
    return SafePreviewReceiptV1(
        **body,
        receipt_sha256=canonical_sha256(body),
        critical_findings=(_FINDING,),
    )


def test_the_preview_keeps_the_bodies_it_hashes():
    source = inspect.getsource(
        __import__(
            "chemsmart.agent.preview", fromlist=["x"]
        ).execute_safe_preview
    )
    assert "critical_findings = tuple(program_validation.findings)" in source


def test_the_validator_carries_the_findings_forward():
    validator = validator_receipt_from_safe_preview(
        node_id="gauss_opt_freq",
        program="gaussian",
        scientific_identity_sha256="8" * 64,
        safe_preview=_failed_preview(),
    )
    assert validator.status == "invalid"
    assert (
        validator.findings
    ), "the digests are meaningless without the bodies they identify"


def test_the_model_is_told_which_field_disagreed():
    validator = validator_receipt_from_safe_preview(
        node_id="gauss_opt_freq",
        program="gaussian",
        scientific_identity_sha256="8" * 64,
        safe_preview=_failed_preview(),
    )
    disclosed = _public_validator_findings(validator)
    rendered = str(disclosed)
    for expected in (
        "gaussian.route.functional_mismatch",
        "route_functional",
        "hf",
        "b3lyp",
    ):
        assert expected in rendered, expected


def test_a_broken_preview_rule_discloses_like_a_validation_finding():
    """A rule objection that renders as blanks is no better than a digest."""

    finding = PreviewRuleFindingV1(
        rule_id="preview.generated_artifact_missing"
    )
    assert finding.field and finding.expected and finding.observed


def test_the_response_the_model_receives_carries_the_findings():
    """Disclosing only to the event log leaves the model where it was."""

    source = inspect.getsource(CommandCompiledToolHostV1._preview_command)
    assert '"critical_findings": _public_validator_findings(validator)' in (
        source
    ), "the model, not only the audit log, has to be able to act on these"


def test_disclosure_does_not_change_receipt_identity():
    """The bodies ride alongside the digest, never inside it.

    ``critical_finding_sha256s`` already binds every body, so integrity is
    unchanged and no receipt recorded before this field loses its identity.
    """

    with_bodies = _failed_preview()
    body = {
        key: value
        for key, value in with_bodies.__dict__.items()
        if key
        not in {
            "receipt_sha256",
            "critical_findings",
            "auxiliary_input_bindings",
        }
    }
    assert with_bodies.receipt_sha256 == canonical_sha256(body)
