from __future__ import annotations

from copy import deepcopy

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
)
from chemsmart.agent.reviews import (
    build_post_grade_review_disposition,
    build_review_disposition,
    build_review_oracle_manifest,
    build_review_oracle_target,
    post_review_experiment_observation,
)
from chemsmart.agent.specialists import CriticReviewRecordV1
from chemsmart.agent.workflows import ScientificReviewFindingV1


_EVIDENCE_A = "a" * 64
_EVIDENCE_B = "b" * 64
_VALIDATOR = "c" * 64
_RAW_GRADE = "d" * 64
_CANDIDATE = {
    "project": {"sp": {"method": "b3lypg", "basis": "def2-svp"}},
    "workflow": {"nodes": ["sp-initial"]},
}


def _finding(
    finding_id: str,
    *,
    rule_id: str = "project.method",
    target_id: str = "candidate.one",
    evidence: tuple[str, ...] = (_EVIDENCE_A,),
) -> ScientificReviewFindingV1:
    body = {
        "schema_version": "chemsmart.scientific-review-finding.v1",
        "finding_id": finding_id,
        "reviewer_role": "critic",
        "rule_id": rule_id,
        "severity": "critical",
        "target_id": target_id,
        "evidence_sha256s": tuple(sorted(set(evidence))),
        "expected": "The frozen host oracle must be satisfied.",
        "observed": "The critic reports a potential mismatch.",
        "disposition": "open",
    }
    return ScientificReviewFindingV1(
        **body, finding_sha256=canonical_sha256(body)
    )


def _review(
    findings: tuple[ScientificReviewFindingV1, ...],
    *,
    candidate: dict | None = None,
    status: str = "complete",
) -> CriticReviewRecordV1:
    candidate = candidate or _CANDIDATE
    body = {
        "schema_version": "chemsmart.critic-review-record.v1",
        "candidate_id": "candidate.one",
        "candidate_sha256": canonical_sha256(candidate),
        "findings": tuple(sorted(findings, key=lambda item: item.finding_id)),
        "status": status,
    }
    return CriticReviewRecordV1(
        **body, review_sha256=canonical_sha256(body)
    )


def _manifest(*targets, coverage=()):
    return build_review_oracle_manifest(
        oracle_id="critic.oracle.one",
        candidate_sha256=canonical_sha256(_CANDIDATE),
        targets=targets,
        validator_coverage_sha256s=coverage,
    )


def _disposition(review, manifest, *, before=None, after=None, verdict="pass"):
    return build_review_disposition(
        candidate_record_before_review=before or _CANDIDATE,
        candidate_record_after_review=after or _CANDIDATE,
        review=review,
        raw_coordinator_grade_sha256=_RAW_GRADE,
        raw_coordinator_verdict=verdict,
        oracle_manifest=manifest,
    )


def test_exact_detection_is_one_true_positive_and_reduces_false_pass():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )

    result = _disposition(
        _review((_finding("finding.one"),)), _manifest(target)
    )

    assert (result.true_positives, result.false_negatives) == (1, 0)
    assert result.overall_recall_basis_points == 10_000
    assert result.critical_recall_basis_points == 10_000
    assert result.raw_false_pass is True
    assert result.review_adjusted_false_pass is False
    assert result.raw_coordinator_verdict == "pass"
    assert result.authority == "post-review-evaluation-only"


def test_missed_defect_remains_a_review_adjusted_false_pass():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )

    result = _disposition(_review(()), _manifest(target))

    assert (result.true_positives, result.false_negatives) == (0, 1)
    assert result.overall_recall_basis_points == 0
    assert result.critical_recall_basis_points == 0
    assert result.raw_false_pass is True
    assert result.review_adjusted_false_pass is True


def test_false_positive_requires_explicit_clean_validator_coverage():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="clean",
        severity="critical",
        registered_evidence_sha256s=(_VALIDATOR,),
        required_validator_sha256s=(_VALIDATOR,),
    )
    with pytest.raises(ContractError, match="incomplete"):
        _manifest(target)

    finding = _finding("finding.clean", evidence=(_VALIDATOR,))
    result = _disposition(
        _review((finding,)), _manifest(target, coverage=(_VALIDATOR,))
    )

    assert (result.false_positives, result.true_negatives) == (1, 0)
    assert result.false_rejection_basis_points == 10_000
    assert result.overall_recall_basis_points is None
    assert result.critical_recall_basis_points is None
    assert result.raw_false_pass is False
    assert result.review_adjusted_false_pass is False


def test_unknown_and_partial_evidence_are_unadjudicated_not_false_positive():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )
    partial = _finding("finding.partial", evidence=(_EVIDENCE_B,))
    unknown = _finding(
        "finding.unknown",
        rule_id="project.unknown",
        evidence=(_EVIDENCE_A,),
    )

    result = _disposition(
        _review((partial, unknown)), _manifest(target)
    )

    assert result.unadjudicated_findings == 2
    assert result.false_positives == 0
    assert result.true_positives == 0
    assert result.false_negatives == 1
    assert {
        item.reason_rule_id for item in result.adjudications
    } == {
        "review.finding.partial-evidence.v1",
        "review.finding.unknown-target.v1",
    }


def test_semantic_duplicate_findings_count_once():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="warning",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )
    findings = (_finding("finding.one"), _finding("finding.two"))

    result = _disposition(_review(findings), _manifest(target))

    assert result.true_positives == 1
    assert result.duplicate_findings == 1
    assert len(result.adjudications) == 1
    assert len(result.adjudications[0].duplicate_finding_sha256s) == 1


def test_digest_tamper_is_rejected_before_scoring():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )
    manifest = _manifest(target)
    object.__setattr__(target, "severity", "info")

    with pytest.raises(ContractError, match="target digest mismatch"):
        _disposition(_review((_finding("finding.one"),)), manifest)


def test_zero_denominator_metrics_are_none_not_zero():
    clean = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="clean",
        severity="warning",
        registered_evidence_sha256s=(_VALIDATOR,),
        required_validator_sha256s=(_VALIDATOR,),
    )
    clean_result = _disposition(
        _review(()), _manifest(clean, coverage=(_VALIDATOR,))
    )
    assert clean_result.overall_recall_basis_points is None
    assert clean_result.critical_recall_basis_points is None
    assert clean_result.false_rejection_basis_points == 0

    warning_defect = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="warning",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )
    defect_result = _disposition(
        _review(()), _manifest(warning_defect)
    )
    assert defect_result.critical_recall_basis_points is None
    assert defect_result.false_rejection_basis_points is None


def test_candidate_byte_identity_and_passive_observation():
    before = deepcopy(_CANDIDATE)
    after = {
        "workflow": {"nodes": ["sp-initial"]},
        "project": {"sp": {"basis": "def2-svp", "method": "b3lypg"}},
    }
    before_bytes = canonical_json(before).encode("utf-8")
    after_bytes = canonical_json(after).encode("utf-8")
    assert before_bytes == after_bytes
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )

    result = _disposition(
        _review((_finding("finding.one"),)),
        _manifest(target),
        before=before,
        after=after,
    )
    observation = post_review_experiment_observation(result)

    assert before == _CANDIDATE
    assert result.candidate_sha256 == canonical_sha256(_CANDIDATE)
    assert observation["authority"] == "post-review-evaluation-only"
    assert "terminal_state" not in canonical_json(observation)
    assert "readiness" not in canonical_json(observation)


def test_candidate_mutation_is_rejected():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )
    changed = deepcopy(_CANDIDATE)
    changed["project"]["sp"]["method"] = "hf"

    with pytest.raises(ContractError, match="changed"):
        _disposition(
            _review((_finding("finding.one"),)),
            _manifest(target),
            before=_CANDIDATE,
            after=changed,
        )


def test_post_grade_analysis_preserves_raw_grade_and_terminal_state():
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="candidate.one",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_EVIDENCE_A,),
    )
    review = _review((_finding("finding.one"),))
    critic = {
        "candidate_sha256": canonical_sha256(_CANDIDATE),
        "candidate_record": canonical_data(_CANDIDATE),
        "review_sha256": review.review_sha256,
        "review": canonical_data(review),
    }
    grade_body = {
        "schema_version": "chemsmart.qwen-pyscf-grade.v1",
        "session_terminal_state": "complete",
        "scientific_state": "previewed",
        "verdict": "pass",
        "checks": (),
    }
    grade = {**grade_body, "grade_sha256": canonical_sha256(grade_body)}
    critic_before = deepcopy(critic)
    grade_before = deepcopy(grade)

    result = build_post_grade_review_disposition(
        critic_observation=critic,
        raw_coordinator_grade_record=grade,
        oracle_manifest=_manifest(target),
    )

    assert critic == critic_before
    assert grade == grade_before
    assert grade["session_terminal_state"] == "complete"
    assert grade["scientific_state"] == "previewed"
    assert grade["verdict"] == "pass"
    assert result.raw_coordinator_grade_sha256 == grade["grade_sha256"]
    assert result.raw_coordinator_verdict == "pass"
    assert result.review_adjusted_false_pass is False
