"""Standalone seeded-defect evaluation for the Qwen/PySCF critic.

The primary D x F x C campaign records whether a fresh critic ran and what it
reported, but it has no answer key for critic efficacy.  This module defines a
separate, bounded study in which candidate bytes, raw coordinator grades, and
defect/clean targets are frozen before a critic session begins.

The critic remains detection-only.  It may inspect registered evidence and
return open findings, but it cannot repair a candidate, approve execution, set
readiness, or replace the raw coordinator grade.  Every score is therefore a
post-review experiment observation rather than a runtime decision.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.reviews import (
    ReviewDispositionV1,
    ReviewOracleManifestV1,
    build_review_disposition,
    build_review_oracle_manifest,
)
from chemsmart.agent.specialists import CriticReviewRecordV1


_STUDY_SCOPE = "standalone-seeded-critic-efficacy"
_PRIMARY_C_SCOPE = "descriptive-only"
_AUTHORITY = "post-review-evaluation-only"
_MAX_CASE_LIMIT = 64
_RAW_VERDICTS = frozenset({"pass", "fail", "inconclusive"})
_READ_ONLY_CRITIC_TOOLS = frozenset(
    {
        "inspect_calculation_artifact",
        "inspect_program_capability",
        "inspect_program_environment",
        "read_project_yaml",
        "validate_project_yaml",
    }
)


@dataclass(frozen=True)
class SeededCriticCaseV1:
    """One preregistered candidate and answer key, frozen before review."""

    schema_version: str
    case_id: str
    benchmark_case_sha256: str
    candidate_sha256: str
    raw_coordinator_grade_sha256: str
    raw_coordinator_verdict: str
    oracle_manifest: ReviewOracleManifestV1
    critic_session_limit: int
    case_contract_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.seeded-critic-case.v1":
            raise ContractError("unsupported seeded critic case schema")
        require_identifier(self.case_id, "case_id")
        for digest, name in (
            (self.benchmark_case_sha256, "benchmark_case_sha256"),
            (self.candidate_sha256, "candidate_sha256"),
            (
                self.raw_coordinator_grade_sha256,
                "raw_coordinator_grade_sha256",
            ),
        ):
            require_sha256(digest, name)
        if self.raw_coordinator_verdict not in _RAW_VERDICTS:
            raise ContractError("unsupported raw coordinator verdict")
        _validate_oracle_manifest(self.oracle_manifest)
        if self.oracle_manifest.candidate_sha256 != self.candidate_sha256:
            raise ContractError("seeded oracle targets another candidate")
        if self.critic_session_limit != 1:
            raise ContractError("seeded critic case requires one fresh critic")
        if self.case_contract_sha256 != canonical_sha256(_case_body(self)):
            raise ContractError("seeded critic case digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(self)


@dataclass(frozen=True)
class SeededCriticStudyManifestV1:
    """Bounded standalone study; it does not upgrade the D x F x C factor."""

    schema_version: str
    study_id: str
    study_scope: str
    primary_dfc_critic_interpretation: str
    case_limit: int
    cases: tuple[SeededCriticCaseV1, ...]
    manifest_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.seeded-critic-study.v1":
            raise ContractError("unsupported seeded critic study schema")
        require_identifier(self.study_id, "study_id")
        if self.study_scope != _STUDY_SCOPE:
            raise ContractError("critic efficacy study is not standalone")
        if self.primary_dfc_critic_interpretation != _PRIMARY_C_SCOPE:
            raise ContractError("primary D x F x C critic must stay descriptive")
        if not 1 <= self.case_limit <= _MAX_CASE_LIMIT:
            raise ContractError("seeded critic case limit is out of range")
        if not self.cases or len(self.cases) > self.case_limit:
            raise ContractError("seeded critic cases exceed the frozen bound")
        order = tuple(item.case_id for item in self.cases)
        if order != tuple(sorted(set(order))):
            raise ContractError("seeded critic cases are not canonical")
        for case in self.cases:
            _validate_case(case)
        if self.manifest_sha256 != canonical_sha256(_manifest_body(self)):
            raise ContractError("seeded critic study digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(self)


@dataclass(frozen=True)
class SeededCriticCaseScoreV1:
    """Detection-only score for exactly one frozen case and critic session."""

    schema_version: str
    authority: str
    case_id: str
    case_contract_sha256: str
    benchmark_case_sha256: str
    review_sha256: str
    critic_tool_names: tuple[str, ...]
    critic_session_count: int
    disposition: ReviewDispositionV1
    score_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.seeded-critic-case-score.v1":
            raise ContractError("unsupported seeded critic score schema")
        if self.authority != _AUTHORITY:
            raise ContractError("critic score cannot own runtime authority")
        require_identifier(self.case_id, "case_id")
        for digest, name in (
            (self.case_contract_sha256, "case_contract_sha256"),
            (self.benchmark_case_sha256, "benchmark_case_sha256"),
            (self.review_sha256, "review_sha256"),
        ):
            require_sha256(digest, name)
        _validate_critic_tools(self.critic_tool_names)
        if self.critic_session_count != 1:
            raise ContractError("critic efficacy requires one fresh session")
        _validate_disposition(self.disposition)
        if self.disposition.review_sha256 != self.review_sha256:
            raise ContractError("critic score review digest mismatch")
        if self.score_sha256 != canonical_sha256(_score_body(self)):
            raise ContractError("seeded critic score digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(self)


@dataclass(frozen=True)
class SeededCriticStudyReportV1:
    """Aggregate efficacy metrics over the exact preregistered case set."""

    schema_version: str
    authority: str
    study_manifest_sha256: str
    primary_dfc_critic_interpretation: str
    case_scores: tuple[SeededCriticCaseScoreV1, ...]
    true_positives: int
    false_positives: int
    false_negatives: int
    true_negatives: int
    critical_true_positives: int
    critical_false_negatives: int
    overall_recall_basis_points: int | None
    critical_recall_basis_points: int | None
    false_rejection_basis_points: int | None
    raw_false_passes: int
    review_adjusted_false_passes: int
    false_pass_reduction_count: int
    false_pass_reduction_basis_points: int | None
    report_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.seeded-critic-report.v1":
            raise ContractError("unsupported seeded critic report schema")
        if self.authority != _AUTHORITY:
            raise ContractError("critic report cannot own runtime authority")
        require_sha256(
            self.study_manifest_sha256, "study_manifest_sha256"
        )
        if self.primary_dfc_critic_interpretation != _PRIMARY_C_SCOPE:
            raise ContractError("primary D x F x C critic must stay descriptive")
        order = tuple(item.case_id for item in self.case_scores)
        if order != tuple(sorted(set(order))) or not order:
            raise ContractError("critic case scores are not canonical")
        for score in self.case_scores:
            _validate_score(score)
        counts = (
            self.true_positives,
            self.false_positives,
            self.false_negatives,
            self.true_negatives,
            self.critical_true_positives,
            self.critical_false_negatives,
            self.raw_false_passes,
            self.review_adjusted_false_passes,
            self.false_pass_reduction_count,
        )
        if any(value < 0 for value in counts):
            raise ContractError("critic efficacy counts cannot be negative")
        if self.review_adjusted_false_passes > self.raw_false_passes:
            raise ContractError("review cannot increase the raw false pass count")
        if self.false_pass_reduction_count != (
            self.raw_false_passes - self.review_adjusted_false_passes
        ):
            raise ContractError("false-pass reduction count is inconsistent")
        for value in (
            self.overall_recall_basis_points,
            self.critical_recall_basis_points,
            self.false_rejection_basis_points,
            self.false_pass_reduction_basis_points,
        ):
            if value is not None and not 0 <= value <= 10_000:
                raise ContractError("critic efficacy basis points are invalid")
        if self.report_sha256 != canonical_sha256(_report_body(self)):
            raise ContractError("seeded critic report digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(self)


def build_seeded_critic_case(
    *,
    case_id: str,
    benchmark_case_sha256: str,
    candidate_record: Mapping[str, Any],
    raw_coordinator_grade_record: Mapping[str, Any],
    oracle_manifest: ReviewOracleManifestV1,
) -> SeededCriticCaseV1:
    """Freeze one case before exposing its candidate to a critic."""

    candidate = canonical_data(dict(candidate_record))
    candidate_sha256 = canonical_sha256(candidate)
    _validate_oracle_manifest(oracle_manifest)
    if oracle_manifest.candidate_sha256 != candidate_sha256:
        raise ContractError("seeded oracle targets another candidate")
    benchmark_digest = require_sha256(
        benchmark_case_sha256, "benchmark_case_sha256"
    )
    grade = _validated_raw_grade(
        raw_coordinator_grade_record,
        expected_case_sha256=benchmark_digest,
    )
    body = {
        "schema_version": "chemsmart.seeded-critic-case.v1",
        "case_id": require_identifier(case_id, "case_id"),
        "benchmark_case_sha256": benchmark_digest,
        "candidate_sha256": candidate_sha256,
        "raw_coordinator_grade_sha256": grade["grade_sha256"],
        "raw_coordinator_verdict": grade["verdict"],
        "oracle_manifest": oracle_manifest,
        "critic_session_limit": 1,
    }
    return SeededCriticCaseV1(
        **body, case_contract_sha256=canonical_sha256(body)
    )


def build_seeded_critic_study_manifest(
    *,
    study_id: str,
    cases: Sequence[SeededCriticCaseV1],
    case_limit: int,
) -> SeededCriticStudyManifestV1:
    """Preregister the complete, bounded case set before any review runs."""

    ordered = tuple(sorted(cases, key=lambda item: item.case_id))
    body = {
        "schema_version": "chemsmart.seeded-critic-study.v1",
        "study_id": require_identifier(study_id, "study_id"),
        "study_scope": _STUDY_SCOPE,
        "primary_dfc_critic_interpretation": _PRIMARY_C_SCOPE,
        "case_limit": int(case_limit),
        "cases": ordered,
    }
    return SeededCriticStudyManifestV1(
        **body, manifest_sha256=canonical_sha256(body)
    )


def score_seeded_critic_case(
    *,
    case: SeededCriticCaseV1,
    observed_benchmark_case_sha256: str,
    candidate_record_before_review: Mapping[str, Any],
    candidate_record_after_review: Mapping[str, Any],
    review: CriticReviewRecordV1,
    critic_tool_names: Sequence[str] = (),
    critic_session_count: int = 1,
) -> SeededCriticCaseScoreV1:
    """Score one review after checking case, authority, and finding digests."""

    _validate_case(case)
    observed_case = require_sha256(
        observed_benchmark_case_sha256,
        "observed_benchmark_case_sha256",
    )
    if observed_case != case.benchmark_case_sha256:
        raise ContractError("critic episode targets another benchmark case")
    tools = tuple(sorted(set(str(item) for item in critic_tool_names)))
    if tuple(critic_tool_names) != tools:
        raise ContractError("critic tool observations are not canonical")
    _validate_critic_tools(tools)
    if int(critic_session_count) != case.critic_session_limit:
        raise ContractError("critic session count differs from preregistration")

    before = canonical_data(dict(candidate_record_before_review))
    after = canonical_data(dict(candidate_record_after_review))
    if canonical_sha256(before) != case.candidate_sha256:
        raise ContractError("critic episode candidate digest mismatch")
    disposition = build_review_disposition(
        candidate_record_before_review=before,
        candidate_record_after_review=after,
        review=review,
        raw_coordinator_grade_sha256=(
            case.raw_coordinator_grade_sha256
        ),
        raw_coordinator_verdict=case.raw_coordinator_verdict,
        oracle_manifest=case.oracle_manifest,
    )
    body = {
        "schema_version": "chemsmart.seeded-critic-case-score.v1",
        "authority": _AUTHORITY,
        "case_id": case.case_id,
        "case_contract_sha256": case.case_contract_sha256,
        "benchmark_case_sha256": case.benchmark_case_sha256,
        "review_sha256": review.review_sha256,
        "critic_tool_names": tools,
        "critic_session_count": 1,
        "disposition": disposition,
    }
    return SeededCriticCaseScoreV1(
        **body, score_sha256=canonical_sha256(body)
    )


def build_seeded_critic_study_report(
    *,
    manifest: SeededCriticStudyManifestV1,
    case_scores: Sequence[SeededCriticCaseScoreV1],
) -> SeededCriticStudyReportV1:
    """Aggregate exact preregistered cases without altering any raw grade."""

    _validate_manifest(manifest)
    scores = tuple(sorted(case_scores, key=lambda item: item.case_id))
    expected = {
        item.case_id: item.case_contract_sha256 for item in manifest.cases
    }
    observed = {item.case_id: item.case_contract_sha256 for item in scores}
    if len(observed) != len(scores) or observed != expected:
        raise ContractError("critic scores do not match frozen study cases")
    for score in scores:
        _validate_score(score)

    dispositions = tuple(item.disposition for item in scores)
    true_positives = sum(item.true_positives for item in dispositions)
    false_positives = sum(item.false_positives for item in dispositions)
    false_negatives = sum(item.false_negatives for item in dispositions)
    true_negatives = sum(item.true_negatives for item in dispositions)
    critical_true_positives = sum(
        item.critical_true_positives for item in dispositions
    )
    critical_false_negatives = sum(
        item.critical_false_negatives for item in dispositions
    )
    raw_false_passes = sum(item.raw_false_pass for item in dispositions)
    adjusted_false_passes = sum(
        item.review_adjusted_false_pass for item in dispositions
    )
    reduction = raw_false_passes - adjusted_false_passes
    body = {
        "schema_version": "chemsmart.seeded-critic-report.v1",
        "authority": _AUTHORITY,
        "study_manifest_sha256": manifest.manifest_sha256,
        "primary_dfc_critic_interpretation": _PRIMARY_C_SCOPE,
        "case_scores": scores,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives,
        "true_negatives": true_negatives,
        "critical_true_positives": critical_true_positives,
        "critical_false_negatives": critical_false_negatives,
        "overall_recall_basis_points": _basis_points(
            true_positives, true_positives + false_negatives
        ),
        "critical_recall_basis_points": _basis_points(
            critical_true_positives,
            critical_true_positives + critical_false_negatives,
        ),
        "false_rejection_basis_points": _basis_points(
            false_positives, false_positives + true_negatives
        ),
        "raw_false_passes": raw_false_passes,
        "review_adjusted_false_passes": adjusted_false_passes,
        "false_pass_reduction_count": reduction,
        "false_pass_reduction_basis_points": _basis_points(
            reduction, raw_false_passes
        ),
    }
    return SeededCriticStudyReportV1(
        **body, report_sha256=canonical_sha256(body)
    )


def _validated_raw_grade(
    value: Mapping[str, Any], *, expected_case_sha256: str
) -> dict[str, Any]:
    grade = canonical_data(dict(value))
    if grade.get("schema_version") != "chemsmart.qwen-pyscf-grade.v1":
        raise ContractError("unsupported raw Qwen/PySCF grade schema")
    grade_sha256 = str(grade.get("grade_sha256") or "")
    require_sha256(grade_sha256, "raw_coordinator_grade_sha256")
    grade_body = {
        key: item for key, item in grade.items() if key != "grade_sha256"
    }
    if grade_sha256 != canonical_sha256(grade_body):
        raise ContractError("raw coordinator grade digest mismatch")
    if str(grade.get("case_sha256") or "") != expected_case_sha256:
        raise ContractError("raw grade targets another benchmark case")
    verdict = str(grade.get("verdict") or "")
    if verdict not in _RAW_VERDICTS:
        raise ContractError("unsupported raw coordinator verdict")
    return grade


def _validate_critic_tools(tools: tuple[str, ...]) -> None:
    if tools != tuple(sorted(set(tools))):
        raise ContractError("critic tool observations are not canonical")
    for tool in tools:
        require_identifier(tool, "critic_tool_name")
    forbidden = set(tools).difference(_READ_ONLY_CRITIC_TOOLS)
    if forbidden:
        raise ContractError(
            "critic used repair, approval, execution, or authority tools"
        )


def _validate_oracle_manifest(manifest: ReviewOracleManifestV1) -> None:
    rebuilt = build_review_oracle_manifest(
        oracle_id=manifest.oracle_id,
        candidate_sha256=manifest.candidate_sha256,
        targets=manifest.targets,
        validator_coverage_sha256s=(
            manifest.validator_coverage_sha256s
        ),
    )
    if rebuilt.manifest_sha256 != manifest.manifest_sha256:
        raise ContractError("seeded oracle manifest digest mismatch")


def _validate_disposition(disposition: ReviewDispositionV1) -> None:
    disposition.__post_init__()
    if disposition.authority != _AUTHORITY:
        raise ContractError("review disposition gained runtime authority")


def _validate_case(case: SeededCriticCaseV1) -> None:
    case.__post_init__()


def _validate_manifest(manifest: SeededCriticStudyManifestV1) -> None:
    manifest.__post_init__()


def _validate_score(score: SeededCriticCaseScoreV1) -> None:
    score.__post_init__()


def _case_body(case: SeededCriticCaseV1) -> dict[str, Any]:
    return {
        key: item
        for key, item in case.__dict__.items()
        if key != "case_contract_sha256"
    }


def _manifest_body(manifest: SeededCriticStudyManifestV1) -> dict[str, Any]:
    return {
        key: item
        for key, item in manifest.__dict__.items()
        if key != "manifest_sha256"
    }


def _score_body(score: SeededCriticCaseScoreV1) -> dict[str, Any]:
    return {
        key: item
        for key, item in score.__dict__.items()
        if key != "score_sha256"
    }


def _report_body(report: SeededCriticStudyReportV1) -> dict[str, Any]:
    return {
        key: item
        for key, item in report.__dict__.items()
        if key != "report_sha256"
    }


def _basis_points(numerator: int, denominator: int) -> int | None:
    if denominator == 0:
        return None
    return round(numerator * 10_000 / denominator)


__all__ = [
    "SeededCriticCaseScoreV1",
    "SeededCriticCaseV1",
    "SeededCriticStudyManifestV1",
    "SeededCriticStudyReportV1",
    "build_seeded_critic_case",
    "build_seeded_critic_study_manifest",
    "build_seeded_critic_study_report",
    "score_seeded_critic_case",
]
