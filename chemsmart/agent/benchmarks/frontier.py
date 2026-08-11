"""Evidence-bound contracts for adapting external frontier benchmarks.

The records in this module bind public task inputs separately from host-only
oracles.  They deliberately contain no benchmark answers and grant no engine
execution authority.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)


_SOURCE_ROLES = frozenset(
    {
        "paper",
        "task-prompt",
        "coordinate",
        "source-scorer",
        "ground-truth",
        "methodology-rubric",
    }
)
_EVIDENCE_CLASSES = frozenset({"preprint", "replication-data"})
_VISIBILITY = frozenset({"public-metadata", "model-input", "host-oracle"})
_GRADING_METHODS = frozenset(
    {"external-deterministic-scorer", "weighted-source-rubric", "all-required"}
)
_ADAPTATION_DECISIONS = frozenset(
    {"adopted", "adapted", "rejected", "no-analogue"}
)


def _body(value: Any, digest_field: str) -> dict[str, Any]:
    return {
        key: item
        for key, item in asdict(value).items()
        if key != digest_field
    }


def _canonical_digests(
    values: tuple[str, ...], field_name: str, *, nonempty: bool = False
) -> None:
    if nonempty and not values:
        raise ContractError(f"{field_name} must not be empty")
    if values != tuple(sorted(set(values))):
        raise ContractError(f"{field_name} must be sorted and unique")
    for value in values:
        require_sha256(value, field_name)


@dataclass(frozen=True)
class BenchmarkSourceArtifactV1:
    """One exact external source artifact and its model-visibility class."""

    schema_version: str
    artifact_id: str
    role: str
    evidence_class: str
    persistent_locator: str
    source_revision: str
    license_spdx: str
    media_type: str
    visibility: str
    content_sha256: str
    size_bytes: int
    record_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.benchmark-source-artifact.v1":
            raise ContractError("unsupported benchmark source artifact schema")
        require_identifier(self.artifact_id, "artifact_id")
        if self.role not in _SOURCE_ROLES:
            raise ContractError("unsupported benchmark source role")
        if self.evidence_class not in _EVIDENCE_CLASSES:
            raise ContractError("unsupported benchmark evidence class")
        if not self.persistent_locator.strip() or not self.source_revision.strip():
            raise ContractError("benchmark source requires locator and revision")
        if not self.license_spdx.strip() or not self.media_type.strip():
            raise ContractError("benchmark source requires license and media type")
        if self.visibility not in _VISIBILITY:
            raise ContractError("unsupported benchmark source visibility")
        require_sha256(self.content_sha256, "content_sha256")
        if self.size_bytes < 1:
            raise ContractError("benchmark source artifact must be non-empty")
        if self.record_sha256 != canonical_sha256(_body(self, "record_sha256")):
            raise ContractError("benchmark source artifact digest mismatch")


@dataclass(frozen=True)
class BenchmarkSourceLedgerV1:
    """Canonical paper and replication-data ledger for one benchmark family."""

    schema_version: str
    ledger_id: str
    benchmark_family: str
    source_locators: tuple[str, ...]
    artifacts: tuple[BenchmarkSourceArtifactV1, ...]
    ledger_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.benchmark-source-ledger.v1":
            raise ContractError("unsupported benchmark source ledger schema")
        require_identifier(self.ledger_id, "ledger_id")
        require_identifier(self.benchmark_family, "benchmark_family")
        if (
            not self.source_locators
            or self.source_locators != tuple(sorted(set(self.source_locators)))
            or any(not item.strip() for item in self.source_locators)
        ):
            raise ContractError("source locators must be sorted, unique, and non-empty")
        keys = tuple(item.artifact_id for item in self.artifacts)
        if not keys or keys != tuple(sorted(set(keys))):
            raise ContractError("source artifacts must be sorted and unique")
        if self.ledger_sha256 != canonical_sha256(_body(self, "ledger_sha256")):
            raise ContractError("benchmark source ledger digest mismatch")

    def artifact(self, artifact_id: str) -> BenchmarkSourceArtifactV1:
        matches = tuple(
            item for item in self.artifacts if item.artifact_id == artifact_id
        )
        if len(matches) != 1:
            raise ContractError("benchmark source artifact is not unique")
        return matches[0]


@dataclass(frozen=True)
class BenchmarkMolecularInputV1:
    """One exact coordinate input and its explicit electronic state."""

    ordinal: int
    system_id: str
    geometry_sha256: str
    charge: int
    multiplicity: int

    def __post_init__(self) -> None:
        if self.ordinal < 1:
            raise ContractError("molecular input ordinal must be positive")
        require_identifier(self.system_id, "system_id")
        require_sha256(self.geometry_sha256, "geometry_sha256")
        if self.multiplicity < 1:
            raise ContractError("molecular input multiplicity must be positive")


@dataclass(frozen=True)
class BenchmarkObservableV1:
    """One required reported observable and its canonical unit."""

    observable_id: str
    unit: str

    def __post_init__(self) -> None:
        require_identifier(self.observable_id, "observable_id")
        require_identifier(self.unit, "unit")


@dataclass(frozen=True)
class BenchmarkChallengeFactorsV1:
    """Answer-free scientific factors defining one borrowed challenge."""

    schema_version: str
    challenge_id: str
    source_ledger_sha256: str
    prompt_artifact_sha256: str
    chemistry_domain: str
    task_family: str
    difficulty_level: str
    systems: tuple[BenchmarkMolecularInputV1, ...]
    fixed_settings: tuple[tuple[str, str], ...]
    stages: tuple[str, ...]
    observables: tuple[BenchmarkObservableV1, ...]
    difficulty_factors: tuple[str, ...]
    challenge_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.benchmark-challenge-factors.v1":
            raise ContractError("unsupported benchmark challenge schema")
        require_identifier(self.challenge_id, "challenge_id")
        require_sha256(self.source_ledger_sha256, "source_ledger_sha256")
        require_sha256(self.prompt_artifact_sha256, "prompt_artifact_sha256")
        require_identifier(self.chemistry_domain, "chemistry_domain")
        require_identifier(self.task_family, "task_family")
        require_identifier(self.difficulty_level, "difficulty_level")
        ordinals = tuple(item.ordinal for item in self.systems)
        if not ordinals or ordinals != tuple(range(1, len(self.systems) + 1)):
            raise ContractError("molecular inputs must preserve consecutive source order")
        ids = tuple(item.system_id for item in self.systems)
        if len(ids) != len(set(ids)):
            raise ContractError("molecular system identifiers must be unique")
        setting_names = tuple(item[0] for item in self.fixed_settings)
        if setting_names != tuple(sorted(set(setting_names))) or any(
            not name.strip() or not value.strip()
            for name, value in self.fixed_settings
        ):
            raise ContractError("fixed settings must be sorted, unique, and non-empty")
        if not self.stages or any(not item.strip() for item in self.stages):
            raise ContractError("benchmark challenge requires ordered stages")
        observable_ids = tuple(item.observable_id for item in self.observables)
        if not observable_ids or observable_ids != tuple(
            sorted(set(observable_ids))
        ):
            raise ContractError("observables must be sorted and unique")
        if self.difficulty_factors != tuple(
            sorted(set(self.difficulty_factors))
        ) or any(not item.strip() for item in self.difficulty_factors):
            raise ContractError("difficulty factors must be sorted and unique")
        if self.challenge_sha256 != canonical_sha256(
            _body(self, "challenge_sha256")
        ):
            raise ContractError("benchmark challenge digest mismatch")


@dataclass(frozen=True)
class BenchmarkRubricCriterionV1:
    """One source or ChemSmart criterion without embedded answer values."""

    criterion_id: str
    channel: str
    weight_basis_points: int
    grading_method: str
    source_locator: str
    required: bool
    criterion_sha256: str

    def __post_init__(self) -> None:
        require_identifier(self.criterion_id, "criterion_id")
        require_identifier(self.channel, "channel")
        if not 0 <= self.weight_basis_points <= 10_000:
            raise ContractError("rubric weight must be between 0 and 10000")
        if self.grading_method not in _GRADING_METHODS:
            raise ContractError("unsupported benchmark grading method")
        if not self.source_locator.strip():
            raise ContractError("rubric criterion requires a source locator")
        if self.criterion_sha256 != canonical_sha256(
            _body(self, "criterion_sha256")
        ):
            raise ContractError("benchmark rubric criterion digest mismatch")


@dataclass(frozen=True)
class BenchmarkRubricProfileV1:
    """Dual source-comparable and ChemSmart-strict scoring profile."""

    schema_version: str
    profile_id: str
    challenge_sha256: str
    source_scorer_sha256: str
    source_ground_truth_sha256: str
    source_methodology_sha256: str
    criteria: tuple[BenchmarkRubricCriterionV1, ...]
    original_channel_weights_basis_points: tuple[tuple[str, int], ...]
    original_aggregation_rule: str
    chemsmart_aggregation_rule: str
    unit_conflict_findings: tuple[str, ...]
    profile_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.benchmark-rubric-profile.v1":
            raise ContractError("unsupported benchmark rubric profile schema")
        require_identifier(self.profile_id, "profile_id")
        for field_name, value in (
            ("challenge_sha256", self.challenge_sha256),
            ("source_scorer_sha256", self.source_scorer_sha256),
            ("source_ground_truth_sha256", self.source_ground_truth_sha256),
            ("source_methodology_sha256", self.source_methodology_sha256),
        ):
            require_sha256(value, field_name)
        ids = tuple(item.criterion_id for item in self.criteria)
        if not ids or ids != tuple(sorted(set(ids))):
            raise ContractError("rubric criteria must be sorted and unique")
        channel_names = tuple(item[0] for item in self.original_channel_weights_basis_points)
        if channel_names != tuple(sorted(set(channel_names))):
            raise ContractError("original rubric channels must be sorted and unique")
        if sum(item[1] for item in self.original_channel_weights_basis_points) != 10_000:
            raise ContractError("original rubric channel weights must sum to 10000")
        for channel, _weight in self.original_channel_weights_basis_points:
            channel_total = sum(
                item.weight_basis_points
                for item in self.criteria
                if item.channel == channel
            )
            if channel_total != 10_000:
                raise ContractError("each original rubric channel must sum to 10000")
        require_identifier(self.original_aggregation_rule, "original_aggregation_rule")
        require_identifier(self.chemsmart_aggregation_rule, "chemsmart_aggregation_rule")
        if self.unit_conflict_findings != tuple(
            sorted(set(self.unit_conflict_findings))
        ):
            raise ContractError("unit conflict findings must be sorted and unique")
        if self.profile_sha256 != canonical_sha256(_body(self, "profile_sha256")):
            raise ContractError("benchmark rubric profile digest mismatch")


@dataclass(frozen=True)
class BenchmarkAdaptationDecisionV1:
    """One explicit adoption, adaptation, rejection, or missing analogue."""

    requirement_id: str
    source_locator: str
    decision: str
    chemsmart_analogue: str
    rationale_rule_id: str
    decision_sha256: str

    def __post_init__(self) -> None:
        require_identifier(self.requirement_id, "requirement_id")
        if not self.source_locator.strip():
            raise ContractError("adaptation decision requires a source locator")
        if self.decision not in _ADAPTATION_DECISIONS:
            raise ContractError("unsupported benchmark adaptation decision")
        if not self.chemsmart_analogue.strip():
            raise ContractError("adaptation decision requires a ChemSmart analogue")
        require_identifier(self.rationale_rule_id, "rationale_rule_id")
        if self.decision_sha256 != canonical_sha256(
            _body(self, "decision_sha256")
        ):
            raise ContractError("benchmark adaptation decision digest mismatch")


@dataclass(frozen=True)
class BenchmarkAdaptationRecordV1:
    """Frozen mapping from the source challenge to ChemSmart semantics."""

    schema_version: str
    adaptation_id: str
    challenge_sha256: str
    rubric_profile_sha256: str
    decisions: tuple[BenchmarkAdaptationDecisionV1, ...]
    execution_scope: str
    adaptation_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.benchmark-adaptation-record.v1":
            raise ContractError("unsupported benchmark adaptation schema")
        require_identifier(self.adaptation_id, "adaptation_id")
        require_sha256(self.challenge_sha256, "challenge_sha256")
        require_sha256(self.rubric_profile_sha256, "rubric_profile_sha256")
        ids = tuple(item.requirement_id for item in self.decisions)
        if not ids or ids != tuple(sorted(set(ids))):
            raise ContractError("adaptation decisions must be sorted and unique")
        if self.execution_scope != "preregistered-not-authorized":
            raise ContractError("frontier case execution scope is not authorized")
        if self.adaptation_sha256 != canonical_sha256(
            _body(self, "adaptation_sha256")
        ):
            raise ContractError("benchmark adaptation digest mismatch")


@dataclass(frozen=True)
class FrontierBenchmarkPreregistrationV1:
    """Sealed binding of public inputs, host oracles, and adaptations."""

    schema_version: str
    preregistration_id: str
    case_id: str
    source_ledger_sha256: str
    challenge_sha256: str
    rubric_profile_sha256: str
    adaptation_sha256: str
    model_input_artifact_sha256s: tuple[str, ...]
    host_oracle_artifact_sha256s: tuple[str, ...]
    held_out: bool
    opening_state: str
    execution_authorized: bool
    preregistration_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.frontier-benchmark-preregistration.v1":
            raise ContractError("unsupported frontier preregistration schema")
        require_identifier(self.preregistration_id, "preregistration_id")
        require_identifier(self.case_id, "case_id")
        for field_name, value in (
            ("source_ledger_sha256", self.source_ledger_sha256),
            ("challenge_sha256", self.challenge_sha256),
            ("rubric_profile_sha256", self.rubric_profile_sha256),
            ("adaptation_sha256", self.adaptation_sha256),
        ):
            require_sha256(value, field_name)
        _canonical_digests(
            self.model_input_artifact_sha256s,
            "model input artifact",
            nonempty=True,
        )
        _canonical_digests(
            self.host_oracle_artifact_sha256s,
            "host oracle artifact",
            nonempty=True,
        )
        if set(self.model_input_artifact_sha256s).intersection(
            self.host_oracle_artifact_sha256s
        ):
            raise ContractError("model inputs and host oracles must be disjoint")
        if not self.held_out or self.opening_state != "sealed":
            raise ContractError("frontier preregistration must remain sealed")
        if self.execution_authorized:
            raise ContractError("frontier preregistration cannot authorize execution")
        if self.preregistration_sha256 != canonical_sha256(
            _body(self, "preregistration_sha256")
        ):
            raise ContractError("frontier preregistration digest mismatch")


def validate_frontier_preregistration(
    *,
    source: BenchmarkSourceLedgerV1,
    challenge: BenchmarkChallengeFactorsV1,
    rubric: BenchmarkRubricProfileV1,
    adaptation: BenchmarkAdaptationRecordV1,
    preregistration: FrontierBenchmarkPreregistrationV1,
) -> None:
    """Cross-check one preregistration without opening host-only answers."""

    if challenge.source_ledger_sha256 != source.ledger_sha256:
        raise ContractError("challenge is bound to another source ledger")
    prompt = tuple(
        item
        for item in source.artifacts
        if item.role == "task-prompt"
        and item.content_sha256 == challenge.prompt_artifact_sha256
    )
    if len(prompt) != 1 or prompt[0].visibility != "model-input":
        raise ContractError("challenge prompt is not one public model input")
    geometry_sha256s = tuple(item.geometry_sha256 for item in challenge.systems)
    coordinate_sha256s = tuple(
        item.content_sha256
        for item in source.artifacts
        if item.role == "coordinate" and item.visibility == "model-input"
    )
    if sorted(geometry_sha256s) != sorted(coordinate_sha256s):
        raise ContractError("challenge geometries differ from the source ledger")
    if rubric.challenge_sha256 != challenge.challenge_sha256:
        raise ContractError("rubric is bound to another challenge")
    role_digests = {
        item.role: item.content_sha256
        for item in source.artifacts
        if item.visibility == "host-oracle"
    }
    if (
        rubric.source_scorer_sha256 != role_digests.get("source-scorer")
        or rubric.source_ground_truth_sha256
        != role_digests.get("ground-truth")
        or rubric.source_methodology_sha256
        != role_digests.get("methodology-rubric")
    ):
        raise ContractError("rubric host oracles differ from the source ledger")
    if (
        adaptation.challenge_sha256 != challenge.challenge_sha256
        or adaptation.rubric_profile_sha256 != rubric.profile_sha256
    ):
        raise ContractError("adaptation is bound to another challenge or rubric")
    if (
        preregistration.source_ledger_sha256 != source.ledger_sha256
        or preregistration.challenge_sha256 != challenge.challenge_sha256
        or preregistration.rubric_profile_sha256 != rubric.profile_sha256
        or preregistration.adaptation_sha256 != adaptation.adaptation_sha256
    ):
        raise ContractError("preregistration ancestry mismatch")
    expected_model_inputs = tuple(
        sorted(
            item.content_sha256
            for item in source.artifacts
            if item.visibility == "model-input"
        )
    )
    expected_host_oracles = tuple(
        sorted(
            item.content_sha256
            for item in source.artifacts
            if item.visibility == "host-oracle"
        )
    )
    if preregistration.model_input_artifact_sha256s != expected_model_inputs:
        raise ContractError("preregistration model-input set mismatch")
    if preregistration.host_oracle_artifact_sha256s != expected_host_oracles:
        raise ContractError("preregistration host-oracle set mismatch")


__all__ = [
    "BenchmarkAdaptationDecisionV1",
    "BenchmarkAdaptationRecordV1",
    "BenchmarkChallengeFactorsV1",
    "BenchmarkMolecularInputV1",
    "BenchmarkObservableV1",
    "BenchmarkRubricCriterionV1",
    "BenchmarkRubricProfileV1",
    "BenchmarkSourceArtifactV1",
    "BenchmarkSourceLedgerV1",
    "FrontierBenchmarkPreregistrationV1",
    "validate_frontier_preregistration",
]
