"""Deterministic, host-owned evaluation of read-only critic findings.

The critic remains an observation-only participant.  This module compares its
immutable review with a preregistered oracle after the coordinator has
finished.  It never changes the candidate, the coordinator grade, scientific
readiness, execution state, or a terminal state.

Adjudication is intentionally narrow: a finding is scoreable only when its
``rule_id``, ``target_id``, and complete evidence digest set exactly match one
frozen oracle target.  Unknown targets and partial evidence remain
``unadjudicated`` and therefore are not counted as false positives.  Clean
targets are scoreable only when their declared deterministic-validator
coverage is complete in the frozen manifest.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.specialists import CriticReviewRecordV1
from chemsmart.agent.workflows import ScientificReviewFindingV1


_TARGET_STATES = frozenset({"defect", "clean"})
_SEVERITIES = frozenset({"info", "warning", "critical"})
_ADJUDICATIONS = frozenset(
    {"true_positive", "false_positive", "unadjudicated"}
)
_RAW_VERDICTS = frozenset({"pass", "fail", "inconclusive"})


@dataclass(frozen=True)
class ReviewOracleTargetV1:
    """One preregistered defect or validator-proven clean review target."""

    schema_version: str
    rule_id: str
    target_id: str
    expected_state: str
    severity: str
    registered_evidence_sha256s: tuple[str, ...]
    required_validator_sha256s: tuple[str, ...]
    target_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.review-oracle-target.v1":
            raise ContractError("unsupported review oracle target schema")
        require_identifier(self.rule_id, "rule_id")
        require_identifier(self.target_id, "target_id")
        if self.expected_state not in _TARGET_STATES:
            raise ContractError("unsupported review oracle target state")
        if self.severity not in _SEVERITIES:
            raise ContractError("unsupported review oracle severity")
        _require_canonical_digests(
            self.registered_evidence_sha256s,
            "registered evidence",
            nonempty=True,
        )
        _require_canonical_digests(
            self.required_validator_sha256s,
            "required validator",
            nonempty=self.expected_state == "clean",
        )
        if not set(self.required_validator_sha256s).issubset(
            self.registered_evidence_sha256s
        ):
            raise ContractError(
                "required validators must be registered target evidence"
            )
        if self.target_sha256 != canonical_sha256(_target_body(self)):
            raise ContractError("review oracle target digest mismatch")


@dataclass(frozen=True)
class ReviewOracleManifestV1:
    """Frozen host oracle for one immutable coordinator candidate."""

    schema_version: str
    oracle_id: str
    candidate_sha256: str
    targets: tuple[ReviewOracleTargetV1, ...]
    validator_coverage_sha256s: tuple[str, ...]
    manifest_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.review-oracle-manifest.v1":
            raise ContractError("unsupported review oracle manifest schema")
        require_identifier(self.oracle_id, "oracle_id")
        require_sha256(self.candidate_sha256, "candidate_sha256")
        _require_canonical_digests(
            self.validator_coverage_sha256s,
            "validator coverage",
        )
        order = tuple(
            (item.rule_id, item.target_id, item.target_sha256)
            for item in self.targets
        )
        if order != tuple(sorted(order)):
            raise ContractError("review oracle targets are not canonical")
        keys = tuple((item.rule_id, item.target_id) for item in self.targets)
        if len(keys) != len(set(keys)):
            raise ContractError("review oracle repeats a rule and target")
        coverage = set(self.validator_coverage_sha256s)
        for target in self.targets:
            _validate_target(target)
            if not set(target.required_validator_sha256s).issubset(coverage):
                raise ContractError(
                    "review oracle has incomplete deterministic-validator "
                    "coverage"
                )
        if self.manifest_sha256 != canonical_sha256(_manifest_body(self)):
            raise ContractError("review oracle manifest digest mismatch")


@dataclass(frozen=True)
class ReviewFindingAdjudicationV1:
    """One deduplicated, post-review finding classification."""

    schema_version: str
    representative_finding_sha256: str
    duplicate_finding_sha256s: tuple[str, ...]
    rule_id: str
    target_id: str
    evidence_sha256s: tuple[str, ...]
    oracle_target_sha256: str
    disposition: str
    reason_rule_id: str
    adjudication_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.review-finding-adjudication.v1":
            raise ContractError("unsupported review adjudication schema")
        require_sha256(
            self.representative_finding_sha256,
            "representative_finding_sha256",
        )
        _require_canonical_digests(
            self.duplicate_finding_sha256s,
            "duplicate finding",
        )
        if self.representative_finding_sha256 in (
            self.duplicate_finding_sha256s
        ):
            raise ContractError("representative finding cannot be a duplicate")
        require_identifier(self.rule_id, "rule_id")
        require_identifier(self.target_id, "target_id")
        _require_canonical_digests(
            self.evidence_sha256s,
            "finding evidence",
        )
        if self.disposition not in _ADJUDICATIONS:
            raise ContractError("unsupported review finding adjudication")
        require_identifier(self.reason_rule_id, "reason_rule_id")
        if self.disposition == "unadjudicated":
            if self.oracle_target_sha256:
                require_sha256(
                    self.oracle_target_sha256, "oracle_target_sha256"
                )
        else:
            require_sha256(
                self.oracle_target_sha256, "oracle_target_sha256"
            )
        if self.adjudication_sha256 != canonical_sha256(
            _adjudication_body(self)
        ):
            raise ContractError("review finding adjudication digest mismatch")


@dataclass(frozen=True)
class ReviewDispositionV1:
    """Analysis-only critic metrics; never a corrected scientific grade."""

    schema_version: str
    authority: str
    candidate_sha256: str
    review_sha256: str
    raw_coordinator_grade_sha256: str
    oracle_manifest_sha256: str
    review_status: str
    raw_coordinator_verdict: str
    adjudications: tuple[ReviewFindingAdjudicationV1, ...]
    missed_defect_target_sha256s: tuple[str, ...]
    unflagged_clean_target_sha256s: tuple[str, ...]
    true_positives: int
    false_positives: int
    false_negatives: int
    true_negatives: int
    critical_true_positives: int
    critical_false_negatives: int
    unadjudicated_findings: int
    duplicate_findings: int
    overall_recall_basis_points: int | None
    critical_recall_basis_points: int | None
    false_rejection_basis_points: int | None
    raw_false_pass: bool
    review_adjusted_false_pass: bool
    disposition_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.review-disposition.v1":
            raise ContractError("unsupported review disposition schema")
        if self.authority != "post-review-evaluation-only":
            raise ContractError(
                "review disposition cannot own runtime authority"
            )
        for digest, name in (
            (self.candidate_sha256, "candidate_sha256"),
            (self.review_sha256, "review_sha256"),
            (
                self.raw_coordinator_grade_sha256,
                "raw_coordinator_grade_sha256",
            ),
            (self.oracle_manifest_sha256, "oracle_manifest_sha256"),
        ):
            require_sha256(digest, name)
        if self.review_status not in {"complete", "blocked", "failed"}:
            raise ContractError("unsupported critic review status")
        if self.raw_coordinator_verdict not in _RAW_VERDICTS:
            raise ContractError("unsupported raw coordinator verdict")
        order = tuple(
            (
                item.rule_id,
                item.target_id,
                item.evidence_sha256s,
                item.representative_finding_sha256,
            )
            for item in self.adjudications
        )
        if order != tuple(sorted(order)):
            raise ContractError("review adjudications are not canonical")
        for item in self.adjudications:
            if item.adjudication_sha256 != canonical_sha256(
                _adjudication_body(item)
            ):
                raise ContractError(
                    "review finding adjudication digest mismatch"
                )
        for values, name in (
            (self.missed_defect_target_sha256s, "missed defect target"),
            (self.unflagged_clean_target_sha256s, "unflagged clean target"),
        ):
            _require_canonical_digests(values, name)
        counts = (
            self.true_positives,
            self.false_positives,
            self.false_negatives,
            self.true_negatives,
            self.critical_true_positives,
            self.critical_false_negatives,
            self.unadjudicated_findings,
            self.duplicate_findings,
        )
        if any(value < 0 for value in counts):
            raise ContractError("review metric counts cannot be negative")
        if self.false_negatives != len(self.missed_defect_target_sha256s):
            raise ContractError("false-negative count differs from targets")
        if self.true_negatives != len(self.unflagged_clean_target_sha256s):
            raise ContractError("true-negative count differs from targets")
        for value in (
            self.overall_recall_basis_points,
            self.critical_recall_basis_points,
            self.false_rejection_basis_points,
        ):
            if value is not None and not 0 <= value <= 10_000:
                raise ContractError("review basis points are out of range")
        if self.review_adjusted_false_pass and not self.raw_false_pass:
            raise ContractError(
                "post-review false pass cannot exceed the raw false pass"
            )
        if self.disposition_sha256 != canonical_sha256(
            _disposition_body(self)
        ):
            raise ContractError("review disposition digest mismatch")

    def public_record(self) -> dict[str, Any]:
        """Return a detached observation suitable for experiment sidecars."""

        return canonical_data(self)


def build_review_oracle_target(
    *,
    rule_id: str,
    target_id: str,
    expected_state: str,
    severity: str,
    registered_evidence_sha256s: Sequence[str],
    required_validator_sha256s: Sequence[str] = (),
) -> ReviewOracleTargetV1:
    """Build one canonical target without consulting critic output."""

    body = {
        "schema_version": "chemsmart.review-oracle-target.v1",
        "rule_id": require_identifier(rule_id, "rule_id"),
        "target_id": require_identifier(target_id, "target_id"),
        "expected_state": str(expected_state),
        "severity": str(severity),
        "registered_evidence_sha256s": tuple(
            sorted(set(registered_evidence_sha256s))
        ),
        "required_validator_sha256s": tuple(
            sorted(set(required_validator_sha256s))
        ),
    }
    return ReviewOracleTargetV1(
        **body, target_sha256=canonical_sha256(body)
    )


def build_review_oracle_manifest(
    *,
    oracle_id: str,
    candidate_sha256: str,
    targets: Sequence[ReviewOracleTargetV1],
    validator_coverage_sha256s: Sequence[str] = (),
) -> ReviewOracleManifestV1:
    """Freeze targets before observing the critic review."""

    ordered = tuple(
        sorted(
            targets,
            key=lambda item: (
                item.rule_id,
                item.target_id,
                item.target_sha256,
            ),
        )
    )
    body = {
        "schema_version": "chemsmart.review-oracle-manifest.v1",
        "oracle_id": require_identifier(oracle_id, "oracle_id"),
        "candidate_sha256": require_sha256(
            candidate_sha256, "candidate_sha256"
        ),
        "targets": ordered,
        "validator_coverage_sha256s": tuple(
            sorted(set(validator_coverage_sha256s))
        ),
    }
    return ReviewOracleManifestV1(
        **body, manifest_sha256=canonical_sha256(body)
    )


def build_review_disposition(
    *,
    candidate_record_before_review: Mapping[str, Any],
    candidate_record_after_review: Mapping[str, Any],
    review: CriticReviewRecordV1,
    raw_coordinator_grade_sha256: str,
    raw_coordinator_verdict: str,
    oracle_manifest: ReviewOracleManifestV1,
) -> ReviewDispositionV1:
    """Score a completed read-only review without changing its inputs.

    The two candidate records are compared as canonical JSON bytes.  The
    returned disposition is an experiment observation only; callers must not
    use it to mutate the raw coordinator grade or any runtime state.
    """

    before_bytes = canonical_json(candidate_record_before_review).encode(
        "utf-8"
    )
    after_bytes = canonical_json(candidate_record_after_review).encode("utf-8")
    if before_bytes != after_bytes:
        raise ContractError("critic changed the canonical candidate bytes")
    candidate_sha256 = canonical_sha256(candidate_record_before_review)
    _validate_manifest(oracle_manifest)
    _validate_review(review)
    if oracle_manifest.candidate_sha256 != candidate_sha256:
        raise ContractError("review oracle targets another candidate")
    if review.candidate_sha256 != candidate_sha256:
        raise ContractError("critic review targets another candidate")
    require_sha256(
        raw_coordinator_grade_sha256,
        "raw_coordinator_grade_sha256",
    )
    if raw_coordinator_verdict not in _RAW_VERDICTS:
        raise ContractError("unsupported raw coordinator verdict")

    target_by_key = {
        (target.rule_id, target.target_id): target
        for target in oracle_manifest.targets
    }
    groups: dict[
        tuple[str, str, tuple[str, ...]],
        list[ScientificReviewFindingV1],
    ] = {}
    for finding in review.findings:
        key = (
            finding.rule_id,
            finding.target_id,
            finding.evidence_sha256s,
        )
        groups.setdefault(key, []).append(finding)

    adjudications: list[ReviewFindingAdjudicationV1] = []
    matched_defects: set[str] = set()
    matched_clean: set[str] = set()
    for key in sorted(groups):
        findings = sorted(
            groups[key], key=lambda item: item.finding_sha256
        )
        representative = findings[0]
        duplicates = tuple(
            item.finding_sha256 for item in findings[1:]
        )
        target = target_by_key.get((key[0], key[1]))
        if target is None:
            disposition = "unadjudicated"
            oracle_target_sha256 = ""
            reason_rule_id = "review.finding.unknown-target.v1"
        elif key[2] != target.registered_evidence_sha256s:
            disposition = "unadjudicated"
            oracle_target_sha256 = target.target_sha256
            reason_rule_id = "review.finding.partial-evidence.v1"
        elif target.expected_state == "defect":
            disposition = "true_positive"
            oracle_target_sha256 = target.target_sha256
            reason_rule_id = "review.finding.true-positive.v1"
            matched_defects.add(target.target_sha256)
        else:
            disposition = "false_positive"
            oracle_target_sha256 = target.target_sha256
            reason_rule_id = "review.finding.false-positive.v1"
            matched_clean.add(target.target_sha256)
        body = {
            "schema_version": (
                "chemsmart.review-finding-adjudication.v1"
            ),
            "representative_finding_sha256": (
                representative.finding_sha256
            ),
            "duplicate_finding_sha256s": duplicates,
            "rule_id": representative.rule_id,
            "target_id": representative.target_id,
            "evidence_sha256s": representative.evidence_sha256s,
            "oracle_target_sha256": oracle_target_sha256,
            "disposition": disposition,
            "reason_rule_id": reason_rule_id,
        }
        adjudications.append(
            ReviewFindingAdjudicationV1(
                **body, adjudication_sha256=canonical_sha256(body)
            )
        )

    defect_targets = tuple(
        item
        for item in oracle_manifest.targets
        if item.expected_state == "defect"
    )
    clean_targets = tuple(
        item
        for item in oracle_manifest.targets
        if item.expected_state == "clean"
    )
    missed_defects = tuple(
        sorted(
            item.target_sha256
            for item in defect_targets
            if item.target_sha256 not in matched_defects
        )
    )
    unflagged_clean = tuple(
        sorted(
            item.target_sha256
            for item in clean_targets
            if item.target_sha256 not in matched_clean
        )
    )
    true_positives = len(matched_defects)
    false_positives = len(matched_clean)
    false_negatives = len(missed_defects)
    true_negatives = len(unflagged_clean)
    critical_defects = tuple(
        item for item in defect_targets if item.severity == "critical"
    )
    critical_true_positives = sum(
        item.target_sha256 in matched_defects for item in critical_defects
    )
    critical_false_negatives = (
        len(critical_defects) - critical_true_positives
    )
    unadjudicated = sum(
        item.disposition == "unadjudicated" for item in adjudications
    )
    duplicates = sum(
        len(item.duplicate_finding_sha256s) for item in adjudications
    )
    raw_false_pass = bool(
        raw_coordinator_verdict == "pass" and defect_targets
    )
    review_adjusted_false_pass = bool(
        raw_coordinator_verdict == "pass" and false_negatives
    )
    body = {
        "schema_version": "chemsmart.review-disposition.v1",
        "authority": "post-review-evaluation-only",
        "candidate_sha256": candidate_sha256,
        "review_sha256": review.review_sha256,
        "raw_coordinator_grade_sha256": raw_coordinator_grade_sha256,
        "oracle_manifest_sha256": oracle_manifest.manifest_sha256,
        "review_status": review.status,
        "raw_coordinator_verdict": raw_coordinator_verdict,
        "adjudications": tuple(adjudications),
        "missed_defect_target_sha256s": missed_defects,
        "unflagged_clean_target_sha256s": unflagged_clean,
        "true_positives": true_positives,
        "false_positives": false_positives,
        "false_negatives": false_negatives,
        "true_negatives": true_negatives,
        "critical_true_positives": critical_true_positives,
        "critical_false_negatives": critical_false_negatives,
        "unadjudicated_findings": unadjudicated,
        "duplicate_findings": duplicates,
        "overall_recall_basis_points": _basis_points(
            true_positives, len(defect_targets)
        ),
        "critical_recall_basis_points": _basis_points(
            critical_true_positives, len(critical_defects)
        ),
        "false_rejection_basis_points": _basis_points(
            false_positives, len(clean_targets)
        ),
        "raw_false_pass": raw_false_pass,
        "review_adjusted_false_pass": review_adjusted_false_pass,
    }
    return ReviewDispositionV1(
        **body, disposition_sha256=canonical_sha256(body)
    )


def post_review_experiment_observation(
    disposition: ReviewDispositionV1,
) -> dict[str, Any]:
    """Wrap a disposition for optional post-review experiment recording.

    This is the intended integration point for ``live_specialists``.  It is
    deliberately not wired into an active campaign automatically: the caller
    must already possess a frozen oracle and raw coordinator grade.
    """

    _validate_disposition(disposition)
    body = {
        "schema_version": "chemsmart.post-review-observation.v1",
        "authority": "post-review-evaluation-only",
        "review_disposition": canonical_data(disposition),
    }
    return {**body, "observation_sha256": canonical_sha256(body)}


def build_post_grade_review_disposition(
    *,
    critic_observation: Mapping[str, Any],
    raw_coordinator_grade_record: Mapping[str, Any],
    oracle_manifest: ReviewOracleManifestV1,
) -> ReviewDispositionV1:
    """Evaluate a persisted critic only after the raw grade is immutable.

    The helper validates and copies the public records, then delegates to the
    same analysis-only disposition builder.  It never returns a replacement
    grade and cannot mutate readiness or terminal state.
    """

    candidate = critic_observation.get("candidate_record")
    review_value = critic_observation.get("review")
    if not isinstance(candidate, Mapping):
        raise ContractError("critic observation lacks its candidate record")
    if not isinstance(review_value, Mapping):
        raise ContractError("critic observation lacks its typed review")
    candidate_record = canonical_data(dict(candidate))
    if str(critic_observation.get("candidate_sha256") or "") != (
        canonical_sha256(candidate_record)
    ):
        raise ContractError("critic observation candidate digest mismatch")
    review = _critic_review_from_record(review_value)
    observed_review_sha256 = str(
        critic_observation.get("review_sha256") or ""
    )
    if observed_review_sha256 != review.review_sha256:
        raise ContractError("critic observation review digest mismatch")

    raw_grade = canonical_data(dict(raw_coordinator_grade_record))
    raw_grade_sha256 = str(raw_grade.get("grade_sha256") or "")
    require_sha256(raw_grade_sha256, "raw_coordinator_grade_sha256")
    raw_grade_body = {
        key: value
        for key, value in raw_grade.items()
        if key != "grade_sha256"
    }
    if raw_grade_sha256 != canonical_sha256(raw_grade_body):
        raise ContractError("raw coordinator grade digest mismatch")
    raw_verdict = str(raw_grade.get("verdict") or "")
    return build_review_disposition(
        candidate_record_before_review=candidate_record,
        candidate_record_after_review=canonical_data(candidate_record),
        review=review,
        raw_coordinator_grade_sha256=raw_grade_sha256,
        raw_coordinator_verdict=raw_verdict,
        oracle_manifest=oracle_manifest,
    )


def _target_body(target: ReviewOracleTargetV1) -> dict[str, Any]:
    return {
        "schema_version": target.schema_version,
        "rule_id": target.rule_id,
        "target_id": target.target_id,
        "expected_state": target.expected_state,
        "severity": target.severity,
        "registered_evidence_sha256s": (
            target.registered_evidence_sha256s
        ),
        "required_validator_sha256s": target.required_validator_sha256s,
    }


def _manifest_body(manifest: ReviewOracleManifestV1) -> dict[str, Any]:
    return {
        "schema_version": manifest.schema_version,
        "oracle_id": manifest.oracle_id,
        "candidate_sha256": manifest.candidate_sha256,
        "targets": manifest.targets,
        "validator_coverage_sha256s": (
            manifest.validator_coverage_sha256s
        ),
    }


def _adjudication_body(
    adjudication: ReviewFindingAdjudicationV1,
) -> dict[str, Any]:
    return {
        "schema_version": adjudication.schema_version,
        "representative_finding_sha256": (
            adjudication.representative_finding_sha256
        ),
        "duplicate_finding_sha256s": (
            adjudication.duplicate_finding_sha256s
        ),
        "rule_id": adjudication.rule_id,
        "target_id": adjudication.target_id,
        "evidence_sha256s": adjudication.evidence_sha256s,
        "oracle_target_sha256": adjudication.oracle_target_sha256,
        "disposition": adjudication.disposition,
        "reason_rule_id": adjudication.reason_rule_id,
    }


def _disposition_body(
    disposition: ReviewDispositionV1,
) -> dict[str, Any]:
    return {
        key: value
        for key, value in disposition.__dict__.items()
        if key != "disposition_sha256"
    }


def _require_canonical_digests(
    values: tuple[str, ...],
    field_name: str,
    *,
    nonempty: bool = False,
) -> None:
    if nonempty and not values:
        raise ContractError(f"{field_name} digests must not be empty")
    if values != tuple(sorted(set(values))):
        raise ContractError(f"{field_name} digests must be sorted and unique")
    for digest in values:
        require_sha256(digest, field_name)


def _validate_target(target: ReviewOracleTargetV1) -> None:
    if target.target_sha256 != canonical_sha256(_target_body(target)):
        raise ContractError("review oracle target digest mismatch")


def _validate_manifest(manifest: ReviewOracleManifestV1) -> None:
    for target in manifest.targets:
        _validate_target(target)
    if manifest.manifest_sha256 != canonical_sha256(_manifest_body(manifest)):
        raise ContractError("review oracle manifest digest mismatch")


def _finding_body(finding: ScientificReviewFindingV1) -> dict[str, Any]:
    return {
        key: value
        for key, value in finding.__dict__.items()
        if key != "finding_sha256"
    }


def _validate_review(review: CriticReviewRecordV1) -> None:
    for finding in review.findings:
        if finding.finding_sha256 != canonical_sha256(_finding_body(finding)):
            raise ContractError("scientific review finding digest mismatch")
    body = {
        "schema_version": review.schema_version,
        "candidate_id": review.candidate_id,
        "candidate_sha256": review.candidate_sha256,
        "findings": review.findings,
        "status": review.status,
    }
    if review.review_sha256 != canonical_sha256(body):
        raise ContractError("critic review digest mismatch")


def _critic_review_from_record(
    value: Mapping[str, Any],
) -> CriticReviewRecordV1:
    raw_findings = value.get("findings")
    if not isinstance(raw_findings, (tuple, list)):
        raise ContractError("critic review record lacks findings")
    findings = []
    for raw in raw_findings:
        if not isinstance(raw, Mapping):
            raise ContractError("critic review finding is malformed")
        fields = {
            key: raw.get(key)
            for key in ScientificReviewFindingV1.__dataclass_fields__
        }
        fields["evidence_sha256s"] = tuple(
            str(item) for item in (raw.get("evidence_sha256s") or ())
        )
        findings.append(ScientificReviewFindingV1(**fields))
    review_fields = {
        "schema_version": str(value.get("schema_version") or ""),
        "candidate_id": str(value.get("candidate_id") or ""),
        "candidate_sha256": str(value.get("candidate_sha256") or ""),
        "findings": tuple(findings),
        "status": str(value.get("status") or ""),
        "review_sha256": str(value.get("review_sha256") or ""),
    }
    return CriticReviewRecordV1(**review_fields)


def _validate_disposition(disposition: ReviewDispositionV1) -> None:
    for item in disposition.adjudications:
        if item.adjudication_sha256 != canonical_sha256(
            _adjudication_body(item)
        ):
            raise ContractError("review finding adjudication digest mismatch")
    if disposition.disposition_sha256 != canonical_sha256(
        _disposition_body(disposition)
    ):
        raise ContractError("review disposition digest mismatch")


def _basis_points(numerator: int, denominator: int) -> int | None:
    if denominator == 0:
        return None
    return round(numerator * 10_000 / denominator)


__all__ = [
    "ReviewDispositionV1",
    "ReviewFindingAdjudicationV1",
    "ReviewOracleManifestV1",
    "ReviewOracleTargetV1",
    "build_post_grade_review_disposition",
    "build_review_disposition",
    "build_review_oracle_manifest",
    "build_review_oracle_target",
    "post_review_experiment_observation",
]
