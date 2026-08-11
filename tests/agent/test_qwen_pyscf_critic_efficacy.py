from __future__ import annotations

from copy import deepcopy

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.experiments.qwen_pyscf_critic_efficacy import (
    build_seeded_critic_case,
    build_seeded_critic_study_manifest,
    build_seeded_critic_study_report,
    score_seeded_critic_case,
)
from chemsmart.agent.reviews import (
    build_review_oracle_manifest,
    build_review_oracle_target,
)
from chemsmart.agent.specialists import CriticReviewRecordV1
from chemsmart.agent.workflows import ScientificReviewFindingV1


_DEFECT_EVIDENCE = "a" * 64
_SECOND_EVIDENCE = "b" * 64
_CLEAN_VALIDATOR = "c" * 64


def _candidate(case_id: str) -> dict:
    return {
        "case_id": case_id,
        "project": {"sp": {"method": "b3lypg", "basis": "def2-svp"}},
        "workflow": {"nodes": ["sp-initial"]},
    }


def _raw_grade(case_sha256: str, *, verdict: str = "pass") -> dict:
    body = {
        "schema_version": "chemsmart.qwen-pyscf-grade.v1",
        "case_sha256": case_sha256,
        "oracle_id": "seeded.critic.oracle",
        "session_terminal_state": "complete",
        "scientific_state": "previewed",
        "verdict": verdict,
        "checks": (),
        "safety_violations": (),
        "successful_tool_calls": 1,
        "failed_tool_calls": 0,
    }
    return {**body, "grade_sha256": canonical_sha256(body)}


def _finding(
    finding_id: str,
    *,
    rule_id: str,
    target_id: str,
    evidence: tuple[str, ...],
) -> ScientificReviewFindingV1:
    body = {
        "schema_version": "chemsmart.scientific-review-finding.v1",
        "finding_id": finding_id,
        "reviewer_role": "critic",
        "rule_id": rule_id,
        "severity": "critical",
        "target_id": target_id,
        "evidence_sha256s": evidence,
        "expected": "The frozen seeded oracle must be satisfied.",
        "observed": "The critic detected a mismatch in the candidate.",
        "disposition": "open",
    }
    return ScientificReviewFindingV1(
        **body, finding_sha256=canonical_sha256(body)
    )


def _review(
    candidate: dict,
    findings: tuple[ScientificReviewFindingV1, ...] = (),
) -> CriticReviewRecordV1:
    body = {
        "schema_version": "chemsmart.critic-review-record.v1",
        "candidate_id": candidate["case_id"],
        "candidate_sha256": canonical_sha256(candidate),
        "findings": tuple(sorted(findings, key=lambda item: item.finding_id)),
        "status": "complete",
    }
    return CriticReviewRecordV1(
        **body, review_sha256=canonical_sha256(body)
    )


def _case(
    case_id: str,
    *,
    expected_state: str = "defect",
    severity: str = "critical",
    evidence: str = _DEFECT_EVIDENCE,
    verdict: str = "pass",
):
    candidate = _candidate(case_id)
    case_sha256 = canonical_sha256({"benchmark_case": case_id})
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id=case_id,
        expected_state=expected_state,
        severity=severity,
        registered_evidence_sha256s=(evidence,),
        required_validator_sha256s=(
            (evidence,) if expected_state == "clean" else ()
        ),
    )
    oracle = build_review_oracle_manifest(
        oracle_id=f"oracle.{case_id}",
        candidate_sha256=canonical_sha256(candidate),
        targets=(target,),
        validator_coverage_sha256s=(
            (evidence,) if expected_state == "clean" else ()
        ),
    )
    case = build_seeded_critic_case(
        case_id=case_id,
        benchmark_case_sha256=case_sha256,
        candidate_record=candidate,
        raw_coordinator_grade_record=_raw_grade(
            case_sha256, verdict=verdict
        ),
        oracle_manifest=oracle,
    )
    return case, candidate


def _score(case, candidate, findings=(), *, tools=()):
    return score_seeded_critic_case(
        case=case,
        observed_benchmark_case_sha256=case.benchmark_case_sha256,
        candidate_record_before_review=candidate,
        candidate_record_after_review=deepcopy(candidate),
        review=_review(candidate, findings),
        critic_tool_names=tools,
    )


def test_seeded_study_reports_detection_and_false_pass_metrics():
    critical, critical_candidate = _case("seeded.critical")
    warning, warning_candidate = _case(
        "seeded.warning", severity="warning", evidence=_SECOND_EVIDENCE
    )
    clean, clean_candidate = _case(
        "seeded.clean",
        expected_state="clean",
        severity="warning",
        evidence=_CLEAN_VALIDATOR,
    )
    critical_finding = _finding(
        "finding.critical",
        rule_id="project.method",
        target_id=critical.case_id,
        evidence=(_DEFECT_EVIDENCE,),
    )
    false_rejection = _finding(
        "finding.clean",
        rule_id="project.method",
        target_id=clean.case_id,
        evidence=(_CLEAN_VALIDATOR,),
    )
    scores = (
        _score(clean, clean_candidate, (false_rejection,)),
        _score(warning, warning_candidate),
        _score(
            critical,
            critical_candidate,
            (critical_finding,),
            tools=("inspect_calculation_artifact",),
        ),
    )
    manifest = build_seeded_critic_study_manifest(
        study_id="qwen.seeded.critic.v1",
        cases=(warning, clean, critical),
        case_limit=3,
    )

    report = build_seeded_critic_study_report(
        manifest=manifest, case_scores=scores
    )

    assert (report.true_positives, report.false_negatives) == (1, 1)
    assert (report.false_positives, report.true_negatives) == (1, 0)
    assert report.overall_recall_basis_points == 5_000
    assert report.critical_recall_basis_points == 10_000
    assert report.false_rejection_basis_points == 10_000
    assert report.raw_false_passes == 2
    assert report.review_adjusted_false_passes == 1
    assert report.false_pass_reduction_count == 1
    assert report.false_pass_reduction_basis_points == 5_000
    assert report.authority == "post-review-evaluation-only"
    assert report.primary_dfc_critic_interpretation == "descriptive-only"


def test_case_and_study_are_frozen_before_review_and_bounded():
    case, _ = _case("seeded.bound")
    manifest = build_seeded_critic_study_manifest(
        study_id="qwen.seeded.bound.v1",
        cases=(case,),
        case_limit=1,
    )

    assert manifest.study_scope == "standalone-seeded-critic-efficacy"
    assert manifest.primary_dfc_critic_interpretation == "descriptive-only"
    assert case.critic_session_limit == 1
    with pytest.raises(ContractError, match="out of range"):
        build_seeded_critic_study_manifest(
            study_id="qwen.seeded.too-small.v1",
            cases=(case,),
            case_limit=0,
        )


def test_case_digest_mismatch_fails_closed_before_scoring():
    case, candidate = _case("seeded.case-digest")

    with pytest.raises(ContractError, match="another benchmark case"):
        score_seeded_critic_case(
            case=case,
            observed_benchmark_case_sha256="f" * 64,
            candidate_record_before_review=candidate,
            candidate_record_after_review=candidate,
            review=_review(candidate),
        )

    object.__setattr__(case, "benchmark_case_sha256", "e" * 64)
    with pytest.raises(ContractError, match="case digest mismatch"):
        score_seeded_critic_case(
            case=case,
            observed_benchmark_case_sha256="e" * 64,
            candidate_record_before_review=candidate,
            candidate_record_after_review=candidate,
            review=_review(candidate),
        )


def test_finding_digest_mismatch_fails_closed_before_classification():
    case, candidate = _case("seeded.finding-digest")
    finding = _finding(
        "finding.tampered",
        rule_id="project.method",
        target_id=case.case_id,
        evidence=(_DEFECT_EVIDENCE,),
    )
    review = _review(candidate, (finding,))
    object.__setattr__(finding, "observed", "Tampered after signing.")

    with pytest.raises(ContractError, match="finding digest mismatch"):
        score_seeded_critic_case(
            case=case,
            observed_benchmark_case_sha256=case.benchmark_case_sha256,
            candidate_record_before_review=candidate,
            candidate_record_after_review=candidate,
            review=review,
        )


@pytest.mark.parametrize(
    ("after_mutation", "tools", "sessions", "message"),
    (
        (True, (), 1, "changed the canonical candidate bytes"),
        (False, ("repair_command",), 1, "authority tools"),
        (False, (), 2, "session count"),
    ),
)
def test_critic_cannot_repair_approve_set_readiness_or_add_sessions(
    after_mutation, tools, sessions, message
):
    case, candidate = _case("seeded.authority")
    after = deepcopy(candidate)
    if after_mutation:
        after["readiness"] = "approved"

    with pytest.raises(ContractError, match=message):
        score_seeded_critic_case(
            case=case,
            observed_benchmark_case_sha256=case.benchmark_case_sha256,
            candidate_record_before_review=candidate,
            candidate_record_after_review=after,
            review=_review(candidate),
            critic_tool_names=tools,
            critic_session_count=sessions,
        )


def test_report_requires_the_exact_preregistered_case_set():
    first, first_candidate = _case("seeded.first")
    second, second_candidate = _case(
        "seeded.second", evidence=_SECOND_EVIDENCE
    )
    manifest = build_seeded_critic_study_manifest(
        study_id="qwen.seeded.exact-set.v1",
        cases=(first, second),
        case_limit=2,
    )
    first_score = _score(first, first_candidate)
    second_score = _score(second, second_candidate)

    with pytest.raises(ContractError, match="frozen study cases"):
        build_seeded_critic_study_report(
            manifest=manifest, case_scores=(first_score,)
        )
    with pytest.raises(ContractError, match="frozen study cases"):
        build_seeded_critic_study_report(
            manifest=manifest,
            case_scores=(first_score, first_score, second_score),
        )


def test_raw_grade_case_and_digest_are_bound_at_preregistration():
    candidate = _candidate("seeded.raw-grade")
    expected_case = canonical_sha256({"case": "expected"})
    wrong_case = canonical_sha256({"case": "wrong"})
    target = build_review_oracle_target(
        rule_id="project.method",
        target_id="seeded.raw-grade",
        expected_state="defect",
        severity="critical",
        registered_evidence_sha256s=(_DEFECT_EVIDENCE,),
    )
    oracle = build_review_oracle_manifest(
        oracle_id="oracle.seeded.raw-grade",
        candidate_sha256=canonical_sha256(candidate),
        targets=(target,),
    )

    with pytest.raises(ContractError, match="another benchmark case"):
        build_seeded_critic_case(
            case_id="seeded.raw-grade",
            benchmark_case_sha256=expected_case,
            candidate_record=candidate,
            raw_coordinator_grade_record=_raw_grade(wrong_case),
            oracle_manifest=oracle,
        )

    grade = _raw_grade(expected_case)
    grade["verdict"] = "fail"
    with pytest.raises(ContractError, match="grade digest mismatch"):
        build_seeded_critic_case(
            case_id="seeded.raw-grade",
            benchmark_case_sha256=expected_case,
            candidate_record=candidate,
            raw_coordinator_grade_record=grade,
            oracle_manifest=oracle,
        )
