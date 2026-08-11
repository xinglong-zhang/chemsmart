"""Deadline-bound controller for Qwen/PySCF harness experiments.

The controller composes :func:`run_live_agent_session`; it does not implement
another provider loop.  Its responsibilities are campaign timing, immutable
episode identity, coordinate isolation, split discipline, deterministic
grading, and a path-free public ledger.  Provider transcripts and Runtime V2
events remain in each episode workspace's private ``.chemsmart-agent`` tree.

There is deliberately no transport-attempt limit.  A caller may supply any
finite or lazily generated development schedule.  The controller stops
dispatching at the exact campaign deadline and never retries a completed
hypothesis.  A rate-limit resume hook may pause before the next unique
episode; a paired retry must be represented by a new repeat index and hence a
new hypothesis ID.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
import json
import math
from pathlib import Path
import re
import shutil
import time
from typing import TYPE_CHECKING, Any, Callable, Iterable, Mapping, Protocol

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
    file_sha256,
    require_sha256,
)
from chemsmart.agent.adaptive_api_campaign import AdaptiveApiCampaignPolicyV1
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    QwenPyscfCaseSpecV1,
    QwenDfcArmV1,
    QwenExperimentPreparationV1,
    QwenPyscfEpisodePlanV1,
    build_case_spec,
    build_episode_plans_from_preparations,
    build_qwen_dfc_arm,
)
from chemsmart.agent.experiments.qwen_pyscf_grading import (
    QwenPyscfDeterministicGradeV1,
    grade_qwen_pyscf_episode,
)
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_MODEL,
    ALIBABA_TOKEN_PLAN_PROVIDER,
    load_agent_provider_selection,
)
from chemsmart.agent.identity import (
    ApprovedMolecularIdentityV1,
    build_approved_molecular_identity,
)
from chemsmart.agent.workflows import HarnessExperimentConfigV1

if TYPE_CHECKING:
    from chemsmart.agent.live_session import CampaignPreparationHostSnapshotV1


_CAMPAIGN_DURATION = timedelta(hours=24)
# The Goal contract reserves its final hour for deterministic reconciliation,
# focused validation, documentation, and local commits.  A new provider
# episode must fit its entire frozen allowance before that reserve begins.
_POST_EPISODE_RESERVE_SECONDS = 3600
_EPISODE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_SPLITS = frozenset({"development", "transfer"})
_TERMINATIONS = frozenset(
    {
        "not_activated",
        "deadline_reached",
        "split_complete",
        "quota_exhausted",
        "credential_revoked",
        "provider_paused",
        "safety_red_line",
    }
)
_PROVIDER_ERROR_CLASSES = frozenset(
    {
        "credential_invalid",
        "quota_exhausted",
        "rate_limited",
        "timeout",
        "transport",
        "provider_5xx",
        "invalid_json",
        "invalid_response_shape",
        "model_unavailable",
        "http_error",
    }
)
_TRANSIENT_ERROR_CLASSES = frozenset(
    {"rate_limited", "timeout", "transport", "provider_5xx"}
)


class _ExperimentIntegrityError(ContractError):
    """A live result did not use the frozen experiment configuration."""


class _ExperimentFactorRealizationError(ContractError):
    """A frozen D/F/C factor was not realized by the live session."""

    def __init__(self, *finding_ids: str) -> None:
        canonical = tuple(sorted(set(finding_ids)))
        if not canonical:
            raise ContractError("factor realization error requires a finding")
        self.finding_ids = canonical
        super().__init__("; ".join(canonical))


def _utc_text(value: datetime) -> str:
    if value.tzinfo is None or value.utcoffset() is None:
        raise ContractError("campaign timestamps must be timezone-aware")
    normalized = value.astimezone(timezone.utc).replace(microsecond=0)
    return normalized.isoformat().replace("+00:00", "Z")


def _parse_utc(value: str) -> datetime:
    text = str(value).strip()
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ContractError("campaign timestamp is not ISO-8601") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ContractError("campaign timestamps must be timezone-aware")
    return parsed.astimezone(timezone.utc)


@dataclass(frozen=True)
class QwenCampaignWindowV1:
    schema_version: str
    campaign_id: str
    activation_utc: str
    deadline_utc: str
    provider: str
    model: str
    transport_attempt_limit: None
    top_up_allowed: bool
    engine_calls: int
    approval_files: int
    window_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-campaign-window.v1":
            raise ContractError("unsupported Qwen campaign window")
        if not self.campaign_id.strip():
            raise ContractError("campaign ID is required")
        activation = _parse_utc(self.activation_utc)
        deadline = _parse_utc(self.deadline_utc)
        if deadline - activation != _CAMPAIGN_DURATION:
            raise ContractError("Qwen campaign window must be exactly 24 hours")
        if (
            self.provider != ALIBABA_TOKEN_PLAN_PROVIDER
            or self.model != ALIBABA_TOKEN_PLAN_MODEL
        ):
            raise ContractError("campaign window must remain Qwen 3.8 Max only")
        if self.transport_attempt_limit is not None or self.top_up_allowed:
            raise ContractError("campaign cannot cap calls or authorize top-up")
        if self.engine_calls or self.approval_files:
            raise ContractError("campaign controller is preview-only")
        body = _without_field(self, "window_sha256")
        if self.window_sha256 != canonical_sha256(body):
            raise ContractError("campaign window digest mismatch")


def build_qwen_campaign_window(
    *, campaign_id: str, activation: datetime
) -> QwenCampaignWindowV1:
    activation_text = _utc_text(activation)
    deadline_text = _utc_text(
        _parse_utc(activation_text) + _CAMPAIGN_DURATION
    )
    body = {
        "schema_version": "chemsmart.qwen-campaign-window.v1",
        "campaign_id": str(campaign_id).strip(),
        "activation_utc": activation_text,
        "deadline_utc": deadline_text,
        "provider": ALIBABA_TOKEN_PLAN_PROVIDER,
        "model": ALIBABA_TOKEN_PLAN_MODEL,
        "transport_attempt_limit": None,
        "top_up_allowed": False,
        "engine_calls": 0,
        "approval_files": 0,
    }
    return QwenCampaignWindowV1(
        **body, window_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ApprovedXyzSourceV1:
    """Runtime-only path plus a path-free approved-source identity."""

    path: Path
    expected_sha256: str
    size_bytes: int
    molecular_identity: ApprovedMolecularIdentityV1 | None
    source_binding_sha256: str

    def __post_init__(self) -> None:
        path = Path(self.path)
        require_sha256(self.expected_sha256, "expected_sha256")
        if not path.is_absolute() or not path.is_file() or path.is_symlink():
            raise ContractError("approved XYZ source must be an absolute regular file")
        if path.suffix.lower() != ".xyz":
            raise ContractError("approved coordinate source must be XYZ")
        if self.size_bytes <= 0 or path.stat().st_size != self.size_bytes:
            raise ContractError("approved XYZ source size changed")
        if file_sha256(path) != self.expected_sha256:
            raise ContractError("approved XYZ source hash changed")
        if (
            self.molecular_identity is not None
            and self.molecular_identity.geometry_sha256 != self.expected_sha256
        ):
            raise ContractError(
                "approved molecular identity targets another XYZ source"
            )
        body = {
            "schema_version": "chemsmart.approved-xyz-source.v1",
            "format": "xyz",
            "sha256": self.expected_sha256,
            "size_bytes": self.size_bytes,
            "approval_scope": "user-approved-coordinate-source",
            "molecular_identity_sha256": (
                self.molecular_identity.identity_sha256
                if self.molecular_identity is not None
                else ""
            ),
        }
        if self.source_binding_sha256 != canonical_sha256(body):
            raise ContractError("approved XYZ source binding mismatch")

    def public_record(self) -> dict[str, Any]:
        return {
            "schema_version": "chemsmart.approved-xyz-source.v1",
            "format": "xyz",
            "sha256": self.expected_sha256,
            "size_bytes": self.size_bytes,
            "approval_scope": "user-approved-coordinate-source",
            "molecular_identity": (
                self.molecular_identity.public_record()
                if self.molecular_identity is not None
                else None
            ),
            "source_binding_sha256": self.source_binding_sha256,
        }


def approved_xyz_source(
    path: str | Path,
    *,
    expected_sha256: str,
    molecular_identity: ApprovedMolecularIdentityV1 | None = None,
) -> ApprovedXyzSourceV1:
    source = Path(path).expanduser()
    if not source.is_absolute():
        raise ContractError("approved XYZ source path must be absolute")
    body = {
        "schema_version": "chemsmart.approved-xyz-source.v1",
        "format": "xyz",
        "sha256": require_sha256(expected_sha256, "expected_sha256"),
        "size_bytes": source.stat().st_size if source.is_file() else 0,
        "approval_scope": "user-approved-coordinate-source",
        "molecular_identity_sha256": (
            molecular_identity.identity_sha256
            if molecular_identity is not None
            else ""
        ),
    }
    return ApprovedXyzSourceV1(
        path=source,
        expected_sha256=body["sha256"],
        size_bytes=body["size_bytes"],
        molecular_identity=molecular_identity,
        source_binding_sha256=canonical_sha256(body),
    )


def approved_xyz_source_from_ledger(
    ledger_path: str | Path, *, artifact_id: str
) -> ApprovedXyzSourceV1:
    """Load one exact coordinate and its identity from an approved-source ledger."""

    ledger = _regular_absolute_file(ledger_path, "coordinate source ledger")
    try:
        document = json.loads(ledger.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ContractError("coordinate source ledger is not valid JSON") from exc
    if document.get("schema_version") != "chemsmart.coordinate-source-ledger.v1":
        raise ContractError("unsupported coordinate source ledger schema")
    matches = tuple(
        row
        for row in document.get("records") or ()
        if isinstance(row, Mapping) and row.get("artifact_id") == artifact_id
    )
    if len(matches) != 1:
        raise ContractError("coordinate source ledger artifact must be unique")
    record = dict(matches[0])
    relative_path = Path(str(record.get("path", "")))
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise ContractError("coordinate source ledger path is unsafe")
    source_path = (ledger.parent / relative_path).resolve()
    expected_sha256 = require_sha256(
        str(record.get("sha256", "")), "coordinate sha256"
    )
    approved_names = tuple(record.get("approved_names") or ())
    if not approved_names:
        approved_names = (str(record.get("identity", "")).strip(),)
    identity = build_approved_molecular_identity(
        identity_id=str(record.get("artifact_id", "")),
        approved_names=approved_names,
        geometry_sha256=expected_sha256,
        coordinate_units=str(record.get("units", "")),
        atom_order=tuple(record.get("atom_order") or ()),
        source_locator=(
            f"{str(record.get('source', '')).strip()} | "
            f"{str(record.get('source_locator', '')).strip()}"
        ),
        source_record_sha256=canonical_sha256(record),
    )
    return approved_xyz_source(
        source_path,
        expected_sha256=expected_sha256,
        molecular_identity=identity,
    )


def bind_case_to_approved_xyz(
    case: QwenPyscfCaseSpecV1, source: ApprovedXyzSourceV1
) -> QwenPyscfCaseSpecV1:
    """Bind one answer-free case to the exact approved coordinate bytes."""

    return build_case_spec(
        case_id=case.case_id,
        split=case.split,
        family=case.family,
        task=case.task,
        expected_observation=case.expected_observation,
        deterministic_oracle_id=case.deterministic_oracle_id,
        source_sha256s=(*case.source_sha256s, source.expected_sha256),
    )


def prepare_qwen_campaign_plans(
    *,
    source: ApprovedXyzSourceV1,
    cases: Iterable[QwenPyscfCaseSpecV1],
    arms: Iterable[QwenDfcArmV1],
    repeats: int,
    provider_config_file: str | Path,
    preparation_workspace_root: str | Path,
    campaign_preparation_snapshot: (
        CampaignPreparationHostSnapshotV1 | None
    ) = None,
) -> tuple[
    tuple[QwenPyscfEpisodePlanV1, ...],
    tuple[QwenExperimentPreparationV1, ...],
    CampaignPreparationHostSnapshotV1,
]:
    """Observe one host snapshot, then construct every freezeable plan.

    The returned snapshot is an explicit runtime dependency for the campaign
    controller.  It is never placed in a module cache and is reused across all
    case/arm/repeat preparations and live episodes.
    """

    if repeats < 1:
        raise ContractError("campaign preparation requires a repetition")
    case_rows = tuple(cases)
    arm_rows = tuple(arms)
    if not case_rows or not arm_rows:
        raise ContractError("campaign preparation requires cases and arms")
    if any(source.expected_sha256 not in case.source_sha256s for case in case_rows):
        raise ContractError("prepared cases must bind the approved XYZ")
    config_path = _regular_absolute_file(
        provider_config_file, "provider config"
    )
    root = Path(preparation_workspace_root).expanduser()
    if not root.is_absolute():
        raise ContractError("preparation workspace root must be absolute")
    root.mkdir(parents=True, exist_ok=True)
    if root.is_symlink():
        raise ContractError("preparation workspace root cannot be a symlink")
    from chemsmart.agent.live_session import (
        build_campaign_preparation_host_snapshot,
        probe_live_experiment_preparation,
        validate_campaign_snapshot_binding,
    )

    host_snapshot = campaign_preparation_snapshot
    if host_snapshot is None:
        snapshot_workspace = root / "campaign-host-snapshot"
        snapshot_workspace.mkdir(mode=0o700, exist_ok=True)
        snapshot_coordinate = snapshot_workspace / "approved.xyz"
        _materialize_exact_coordinate(
            source=source,
            workspace=snapshot_workspace,
            coordinate=snapshot_coordinate,
            role="campaign host snapshot",
        )
        host_snapshot = build_campaign_preparation_host_snapshot(
            provider=ALIBABA_TOKEN_PLAN_PROVIDER,
            provider_config_file=config_path,
            workspace=snapshot_workspace,
            approved_molecular_identity=source.molecular_identity,
        )
    else:
        validate_campaign_snapshot_binding(
            snapshot=host_snapshot,
            provider_config_file=config_path,
            provider=ALIBABA_TOKEN_PLAN_PROVIDER,
            artifact_sha256=source.expected_sha256,
            approved_identity_sha256=(
                source.molecular_identity.identity_sha256
                if source.molecular_identity is not None
                else ""
            ),
        )
    preparations = []
    for repeat_index in range(repeats):
        for case in sorted(case_rows, key=lambda item: item.case_id):
            for arm in sorted(arm_rows, key=lambda item: item.arm_id):
                episode_id = f"{case.case_id}.{arm.arm_id}.r{repeat_index}"
                workspace = root / episode_id
                workspace.mkdir(mode=0o700, exist_ok=True)
                coordinate = workspace / "approved.xyz"
                _materialize_exact_coordinate(
                    source=source,
                    workspace=workspace,
                    coordinate=coordinate,
                    role="preparation",
                )
                preparations.append(
                    probe_live_experiment_preparation(
                        task=case.task,
                        provider=ALIBABA_TOKEN_PLAN_PROVIDER,
                        provider_config_file=config_path,
                        workspace=workspace,
                        experiment_arm=arm,
                        experiment_case=case,
                        experiment_repeat_index=repeat_index,
                        approved_molecular_identity=source.molecular_identity,
                        campaign_preparation_snapshot=host_snapshot,
                    )
                )
    preparation_rows = tuple(preparations)
    plans = build_episode_plans_from_preparations(
        cases=case_rows,
        arms=arm_rows,
        preparations=preparation_rows,
    )
    return plans, preparation_rows, host_snapshot


def _materialize_exact_coordinate(
    *,
    source: ApprovedXyzSourceV1,
    workspace: Path,
    coordinate: Path,
    role: str,
) -> None:
    if coordinate.exists():
        if coordinate.is_symlink() or not coordinate.is_file():
            raise ContractError(f"{role} coordinate target is unsafe")
        if file_sha256(coordinate) != source.expected_sha256:
            raise ContractError(f"{role} workspace contains another geometry")
    else:
        shutil.copyfile(source.path, coordinate)
    if tuple(sorted(workspace.glob("*.xyz"))) != (coordinate,):
        raise ContractError(f"{role} workspace must contain exactly one XYZ")


@dataclass(frozen=True)
class FrozenTransferManifestV1:
    schema_version: str
    campaign_window_sha256: str
    approved_source_binding_sha256: str
    development_case_sha256s: tuple[str, ...]
    transfer_case_sha256s: tuple[str, ...]
    development_configuration_sha256s: tuple[str, ...]
    transfer_plan_sha256s: tuple[str, ...]
    manifest_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-transfer-freeze.v1":
            raise ContractError("unsupported Qwen transfer freeze")
        for name in (
            "campaign_window_sha256",
            "approved_source_binding_sha256",
        ):
            require_sha256(getattr(self, name), name)
        for name in (
            "development_case_sha256s",
            "transfer_case_sha256s",
            "development_configuration_sha256s",
            "transfer_plan_sha256s",
        ):
            values = getattr(self, name)
            if not values or values != tuple(sorted(set(values))):
                raise ContractError(f"{name} must be non-empty and canonical")
            for digest in values:
                require_sha256(digest, name)
        if set(self.development_case_sha256s) & set(
            self.transfer_case_sha256s
        ):
            raise ContractError("development and transfer cases overlap")
        body = _without_field(self, "manifest_sha256")
        if self.manifest_sha256 != canonical_sha256(body):
            raise ContractError("Qwen transfer freeze digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(asdict(self))


def build_frozen_transfer_manifest(
    *,
    window: QwenCampaignWindowV1,
    source: ApprovedXyzSourceV1,
    cases: Iterable[QwenPyscfCaseSpecV1],
    plans: Iterable[QwenPyscfEpisodePlanV1],
) -> FrozenTransferManifestV1:
    case_rows = tuple(cases)
    plan_rows = tuple(plans)
    by_sha = {item.case_sha256: item for item in case_rows}
    if len(by_sha) != len(case_rows):
        raise ContractError("campaign cases are not unique")
    if any(source.expected_sha256 not in item.source_sha256s for item in case_rows):
        raise ContractError("every campaign case must bind the approved XYZ")
    unknown = tuple(item for item in plan_rows if item.case_sha256 not in by_sha)
    if unknown:
        raise ContractError("campaign plan references an unknown case")
    development_plans = tuple(
        item for item in plan_rows if by_sha[item.case_sha256].split == "development"
    )
    transfer_plans = tuple(
        item for item in plan_rows if by_sha[item.case_sha256].split == "transfer"
    )
    if not development_plans or not transfer_plans:
        raise ContractError("freeze requires development and transfer plans")
    body = {
        "schema_version": "chemsmart.qwen-transfer-freeze.v1",
        "campaign_window_sha256": window.window_sha256,
        "approved_source_binding_sha256": source.source_binding_sha256,
        "development_case_sha256s": tuple(
            sorted(
                item.case_sha256
                for item in case_rows
                if item.split == "development"
            )
        ),
        "transfer_case_sha256s": tuple(
            sorted(
                item.case_sha256 for item in case_rows if item.split == "transfer"
            )
        ),
        "development_configuration_sha256s": tuple(
            sorted(
                {
                    _experiment_signature(item)
                    for item in development_plans
                }
            )
        ),
        "transfer_plan_sha256s": tuple(
            sorted(item.plan_sha256 for item in transfer_plans)
        ),
    }
    return FrozenTransferManifestV1(
        **body, manifest_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class QwenPyscfEpisodeLedgerV1:
    schema_version: str
    episode_id: str
    split: str
    case_sha256: str
    plan_sha256: str
    hypothesis_sha256: str
    experiment_config_sha256: str
    source_sha256: str
    workspace_binding_sha256: str
    result_sha256: str
    grade_sha256: str
    session_terminal_state: str
    scientific_state: str
    verdict: str
    failure_class: str
    factor_realization_status: str
    factor_realization_findings: tuple[str, ...]
    safety_violations: tuple[str, ...]
    ledger_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-episode-ledger.v1":
            raise ContractError("unsupported Qwen episode ledger")
        if self.split not in _SPLITS:
            raise ContractError("episode ledger split is invalid")
        for name in (
            "case_sha256",
            "plan_sha256",
            "hypothesis_sha256",
            "experiment_config_sha256",
            "source_sha256",
            "workspace_binding_sha256",
            "result_sha256",
            "grade_sha256",
        ):
            require_sha256(getattr(self, name), name)
        if self.verdict not in {"pass", "fail", "inconclusive"}:
            raise ContractError("episode ledger verdict is invalid")
        if self.factor_realization_status not in {
            "valid",
            "invalid",
            "not_observed",
        }:
            raise ContractError("factor realization status is invalid")
        if self.factor_realization_findings != tuple(
            sorted(set(self.factor_realization_findings))
        ):
            raise ContractError("factor realization findings are not canonical")
        if (
            self.factor_realization_status == "invalid"
        ) != bool(self.factor_realization_findings):
            raise ContractError(
                "invalid factor realization must have findings"
            )
        if (
            self.factor_realization_status == "invalid"
        ) != (self.failure_class == "experiment_factor_invalid"):
            raise ContractError(
                "factor realization status and failure class differ"
            )
        if self.safety_violations != tuple(sorted(set(self.safety_violations))):
            raise ContractError("episode safety findings are not canonical")
        body = _without_field(self, "ledger_sha256")
        if self.ledger_sha256 != canonical_sha256(body):
            raise ContractError("episode ledger digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(asdict(self))


@dataclass(frozen=True)
class QwenPyscfSplitLedgerV1:
    schema_version: str
    campaign_window_sha256: str
    freeze_manifest_sha256: str
    split: str
    started_utc: str
    observed_at_utc: str
    termination_reason: str
    episode_ledgers: tuple[QwenPyscfEpisodeLedgerV1, ...]
    last_hypothesis_sha256: str
    ledger_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-split-ledger.v1":
            raise ContractError("unsupported Qwen split ledger")
        require_sha256(self.campaign_window_sha256, "campaign_window_sha256")
        require_sha256(self.freeze_manifest_sha256, "freeze_manifest_sha256")
        if self.split not in _SPLITS:
            raise ContractError("split ledger split is invalid")
        _parse_utc(self.started_utc)
        _parse_utc(self.observed_at_utc)
        if self.termination_reason not in _TERMINATIONS:
            raise ContractError("split ledger termination is invalid")
        ids = tuple(item.episode_id for item in self.episode_ledgers)
        if len(ids) != len(set(ids)):
            raise ContractError("split ledger repeats an episode")
        if any(item.split != self.split for item in self.episode_ledgers):
            raise ContractError("split ledger contains another split")
        if self.last_hypothesis_sha256:
            require_sha256(
                self.last_hypothesis_sha256, "last_hypothesis_sha256"
            )
        body = _without_field(self, "ledger_sha256")
        if self.ledger_sha256 != canonical_sha256(body):
            raise ContractError("split ledger digest mismatch")

    def public_summary_json(self) -> str:
        """Return a path-free, reasoning-free public campaign record."""

        return canonical_json(asdict(self))


class LiveEpisodeRunner(Protocol):
    def __call__(self, **kwargs: Any) -> Any: ...


class EpisodeOutcomeObserver(Protocol):
    def __call__(
        self,
        *,
        plan: QwenPyscfEpisodePlanV1,
        case: QwenPyscfCaseSpecV1,
        result: Any | None,
        grade: QwenPyscfDeterministicGradeV1,
        ledger: QwenPyscfEpisodeLedgerV1,
    ) -> None: ...


class RetryResumeHook(Protocol):
    def __call__(
        self,
        *,
        episode_id: str,
        failure_class: str,
        observed_retry_after_seconds: float | None,
        seconds_until_deadline: float,
    ) -> float | None: ...


def _default_live_runner(**kwargs: Any) -> Any:
    from chemsmart.agent.live_session import run_live_agent_session

    return run_live_agent_session(**kwargs)


def _default_resume_hook(
    *,
    episode_id: str,
    failure_class: str,
    observed_retry_after_seconds: float | None,
    seconds_until_deadline: float,
) -> float | None:
    del episode_id, failure_class
    if observed_retry_after_seconds is None:
        return None
    delay = float(observed_retry_after_seconds)
    return delay if 0 <= delay < seconds_until_deadline else None


@dataclass(frozen=True)
class _PreparedEpisode:
    plan: QwenPyscfEpisodePlanV1
    case: QwenPyscfCaseSpecV1
    workspace: Path
    workspace_binding_sha256: str


@dataclass(frozen=True)
class _EpisodeOutcome:
    ledger: QwenPyscfEpisodeLedgerV1
    resume_after_seconds: float | None


class QwenPyscfCampaignControllerV1:
    """Compose preview-only Qwen episodes under one exact 24-hour window."""

    def __init__(
        self,
        *,
        window: QwenCampaignWindowV1,
        freeze: FrozenTransferManifestV1,
        source: ApprovedXyzSourceV1,
        cases: Iterable[QwenPyscfCaseSpecV1],
        provider_config_file: str | Path,
        secret_file: str | Path,
        workspace_root: str | Path,
        campaign_preparation_snapshot: (
            CampaignPreparationHostSnapshotV1 | None
        ) = None,
        max_concurrency: int = 1,
        runner: LiveEpisodeRunner = _default_live_runner,
        outcome_observer: EpisodeOutcomeObserver | None = None,
        authorized_transfer_plan_sha256s: Iterable[str] = (),
        clock: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
        sleeper: Callable[[float], None] = time.sleep,
        resume_hook: RetryResumeHook = _default_resume_hook,
    ) -> None:
        if freeze.campaign_window_sha256 != window.window_sha256:
            raise ContractError("transfer freeze belongs to another campaign")
        if freeze.approved_source_binding_sha256 != source.source_binding_sha256:
            raise ContractError("transfer freeze binds another coordinate source")
        if not 1 <= int(max_concurrency) <= 4:
            raise ContractError("campaign concurrency must be within [1, 4]")
        self.window = window
        self.freeze = freeze
        self.source = source
        case_rows = tuple(cases)
        self.cases = {item.case_sha256: item for item in case_rows}
        if not case_rows:
            raise ContractError("campaign controller requires cases")
        if len(self.cases) != len(case_rows):
            raise ContractError("campaign cases are not unique")
        expected_cases = set(freeze.development_case_sha256s) | set(
            freeze.transfer_case_sha256s
        )
        if set(self.cases) != expected_cases:
            raise ContractError("controller cases differ from the frozen registry")
        self.provider_config_file = _regular_absolute_file(
            provider_config_file, "provider config"
        )
        self.secret_file = _regular_absolute_file(secret_file, "secret file")
        self.workspace_root = Path(workspace_root).expanduser()
        if not self.workspace_root.is_absolute():
            raise ContractError("campaign workspace root must be absolute")
        self.workspace_root.mkdir(parents=True, exist_ok=True)
        if self.workspace_root.is_symlink():
            raise ContractError("campaign workspace root cannot be a symlink")
        self.max_concurrency = int(max_concurrency)
        self.runner = runner
        self.outcome_observer = outcome_observer
        self.authorized_transfer_plan_sha256s = frozenset(
            require_sha256(item, "authorized transfer plan sha256")
            for item in authorized_transfer_plan_sha256s
        )
        if not self.authorized_transfer_plan_sha256s.issubset(
            freeze.transfer_plan_sha256s
        ):
            raise ContractError("transfer authorization is outside the freeze")
        self.clock = clock
        self.sleeper = sleeper
        self.resume_hook = resume_hook
        self.campaign_preparation_snapshot = campaign_preparation_snapshot
        self.policy = AdaptiveApiCampaignPolicyV1()
        self._dispatched_episode_ids: set[str] = set()
        self._validate_qwen_only_profile()

    def run_split(
        self,
        *,
        split: str,
        plans: Iterable[QwenPyscfEpisodePlanV1],
    ) -> QwenPyscfSplitLedgerV1:
        """Run one split; transfer plans must match the preregistered freeze.

        ``plans`` is intentionally an iterable rather than a count-bounded
        sequence.  A development plan source may therefore yield adaptive,
        uniquely identified episodes until the deadline.  Transfer remains a
        finite frozen set.
        """

        split = str(split).strip().lower()
        if split not in _SPLITS:
            raise ContractError("campaign split must be development or transfer")
        started = self._now()
        activation = _parse_utc(self.window.activation_utc)
        deadline = _parse_utc(self.window.deadline_utc)
        if started < activation:
            return self._split_ledger(
                split=split,
                started=started,
                observed=started,
                termination="not_activated",
                ledgers=(),
            )
        if started >= deadline:
            return self._split_ledger(
                split=split,
                started=started,
                observed=started,
                termination="deadline_reached",
                ledgers=(),
            )

        iterator = iter(plans)
        seen_episode_ids: set[str] = set()
        ledgers: list[QwenPyscfEpisodeLedgerV1] = []
        termination = "split_complete"
        while self._now() < deadline:
            batch: list[_PreparedEpisode] = []
            for _ in range(self.max_concurrency):
                if self._now() >= deadline:
                    termination = "deadline_reached"
                    break
                try:
                    plan = next(iterator)
                except StopIteration:
                    break
                remaining_seconds = (deadline - self._now()).total_seconds()
                if remaining_seconds <= (
                    plan.experiment_config.wall_time_seconds
                    + _POST_EPISODE_RESERVE_SECONDS
                ):
                    # A frozen episode may use its entire wall-time allowance.
                    # Do not begin work that the exact campaign window cannot
                    # contain; the remaining time is reserved for local
                    # reconciliation and final reporting.
                    termination = "deadline_reached"
                    break
                prepared = self._prepare_episode(
                    plan=plan,
                    split=split,
                    seen_episode_ids=seen_episode_ids,
                )
                batch.append(prepared)
            if not batch:
                break
            if self._now() >= deadline:
                termination = "deadline_reached"
                break
            outcomes = self._run_independent_batch(batch)
            ledgers.extend(outcome.ledger for outcome in outcomes)
            requested_pause_seconds: list[float] = []
            batch_termination = ""
            for outcome in outcomes:
                failure = outcome.ledger.failure_class
                if outcome.ledger.safety_violations:
                    batch_termination = "safety_red_line"
                    continue
                if failure == "experiment_integrity_violation":
                    batch_termination = "safety_red_line"
                    continue
                if batch_termination == "safety_red_line":
                    continue
                if failure == "credential_revoked":
                    batch_termination = "credential_revoked"
                    continue
                if batch_termination == "credential_revoked":
                    continue
                if failure in {
                    "quota_exhausted",
                    "rate_limited",
                    "transient_transport",
                }:
                    delay = outcome.resume_after_seconds
                    remaining = (deadline - self._now()).total_seconds()
                    if delay is not None and delay < (
                        remaining - _POST_EPISODE_RESERVE_SECONDS
                    ):
                        requested_pause_seconds.append(delay)
                    elif failure == "quota_exhausted":
                        if not batch_termination:
                            batch_termination = "quota_exhausted"
                    elif not batch_termination:
                        batch_termination = "provider_paused"
            if batch_termination:
                termination = batch_termination
                break
            if requested_pause_seconds:
                self.sleeper(max(requested_pause_seconds))
        observed = self._now()
        if observed >= deadline and termination == "split_complete":
            termination = "deadline_reached"
        return self._split_ledger(
            split=split,
            started=started,
            observed=observed,
            termination=termination,
            ledgers=tuple(ledgers),
        )

    def _validate_qwen_only_profile(self) -> None:
        selection = load_agent_provider_selection(
            self.provider_config_file,
            requested_profile=ALIBABA_TOKEN_PLAN_PROVIDER,
        )
        profile = selection.active_profile
        if (
            profile.provider != ALIBABA_TOKEN_PLAN_PROVIDER
            or profile.model != ALIBABA_TOKEN_PLAN_MODEL
            or selection.fallback_profiles
        ):
            raise ContractError(
                "campaign requires Qwen 3.8 Max with no fallback profiles"
            )
        if self.campaign_preparation_snapshot is not None:
            from chemsmart.agent.live_session import (
                validate_campaign_snapshot_binding,
            )

            validate_campaign_snapshot_binding(
                snapshot=self.campaign_preparation_snapshot,
                provider_config_file=self.provider_config_file,
                provider=ALIBABA_TOKEN_PLAN_PROVIDER,
                artifact_sha256=self.source.expected_sha256,
                approved_identity_sha256=(
                    self.source.molecular_identity.identity_sha256
                    if self.source.molecular_identity is not None
                    else ""
                ),
            )

    def _prepare_episode(
        self,
        *,
        plan: QwenPyscfEpisodePlanV1,
        split: str,
        seen_episode_ids: set[str],
    ) -> _PreparedEpisode:
        if not isinstance(plan, QwenPyscfEpisodePlanV1):
            raise ContractError("campaign plan is not typed")
        if _EPISODE_ID.fullmatch(plan.episode_id) is None:
            raise ContractError("episode ID is not workspace-safe")
        if plan.episode_id in seen_episode_ids:
            raise ContractError("campaign episode identity is not unique")
        case = self.cases.get(plan.case_sha256)
        if case is None or case.split != split:
            raise ContractError("campaign plan crosses the selected split")
        feedback_label = (
            "causal"
            if plan.experiment_config.feedback_projection == "causal-v1"
            else "full"
        )
        expected_arm_id = (
            f"d{int(plan.experiment_config.decomposition)}-"
            f"f{feedback_label}-"
            f"c{int(plan.experiment_config.critic)}"
        )
        expected_episode_id = (
            f"{case.case_id}.{expected_arm_id}.r{plan.repeat_index}"
        )
        if (
            plan.episode_id != expected_episode_id
            or plan.experiment_config.experiment_id != expected_episode_id.lower()
            or plan.hypothesis.changed_factor != "D/F/C=" + expected_arm_id
        ):
            raise ContractError(
                "episode identity must be the unique case/arm/repeat tuple"
            )
        if plan.episode_id in self._dispatched_episode_ids:
            raise ContractError("campaign already dispatched this hypothesis")
        if self.source.expected_sha256 not in plan.hypothesis.source_sha256s:
            raise ContractError("episode hypothesis omits the approved XYZ hash")
        if self.campaign_preparation_snapshot is not None and (
            self.campaign_preparation_snapshot.snapshot_sha256
            not in plan.hypothesis.source_sha256s
        ):
            raise ContractError("episode hypothesis omits the host snapshot")
        if plan.engine_calls != 0:
            raise ContractError("campaign episode cannot call a chemistry engine")
        if plan.experiment_config.provider_id != ALIBABA_TOKEN_PLAN_PROVIDER:
            raise ContractError("campaign episode is not Qwen-only")
        if plan.experiment_config.max_concurrency < self.max_concurrency:
            raise ContractError("episode does not authorize campaign concurrency")
        if split == "transfer":
            if plan.plan_sha256 not in self.freeze.transfer_plan_sha256s:
                raise ContractError("transfer plan was not frozen")
            if plan.plan_sha256 not in self.authorized_transfer_plan_sha256s:
                raise ContractError("transfer plan has not been opened")
        elif (
            _experiment_signature(plan)
            not in self.freeze.development_configuration_sha256s
        ):
            raise ContractError("development configuration was not frozen")
        workspace = self.workspace_root / plan.episode_id
        workspace.mkdir(mode=0o700, parents=False, exist_ok=True)
        if workspace.is_symlink():
            raise ContractError("episode workspace cannot be a symlink")
        coordinate = workspace / "approved.xyz"
        if coordinate.exists():
            if not coordinate.is_file() or coordinate.is_symlink():
                raise ContractError("episode coordinate target is unsafe")
            if file_sha256(coordinate) != self.source.expected_sha256:
                raise ContractError("episode workspace contains another geometry")
        else:
            shutil.copyfile(self.source.path, coordinate)
        if file_sha256(coordinate) != self.source.expected_sha256:
            raise ContractError("copied episode geometry failed its hash check")
        xyz_files = tuple(sorted(workspace.glob("*.xyz")))
        if xyz_files != (coordinate,):
            raise ContractError("episode workspace must contain exactly one XYZ")
        seen_episode_ids.add(plan.episode_id)
        self._dispatched_episode_ids.add(plan.episode_id)
        binding = canonical_sha256(
            {
                "schema_version": "chemsmart.episode-workspace-binding.v1",
                "episode_id": plan.episode_id,
                "source_binding_sha256": self.source.source_binding_sha256,
                "coordinate_sha256": self.source.expected_sha256,
                "mutable_owner": plan.episode_id,
            }
        )
        return _PreparedEpisode(
            plan=plan,
            case=case,
            workspace=workspace,
            workspace_binding_sha256=binding,
        )

    def _run_independent_batch(
        self, batch: list[_PreparedEpisode]
    ) -> tuple[_EpisodeOutcome, ...]:
        owners = tuple(item.plan.episode_id for item in batch)
        if len(owners) != len(set(owners)):
            raise ContractError("concurrent episodes share a mutable owner")
        if len(batch) == 1:
            return (self._run_episode(batch[0]),)
        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = [pool.submit(self._run_episode, item) for item in batch]
            return tuple(future.result() for future in futures)

    def _run_episode(self, item: _PreparedEpisode) -> _EpisodeOutcome:
        plan = item.plan
        arm_config = plan.experiment_config
        arm = build_qwen_dfc_arm(
            decomposition=arm_config.decomposition,
            feedback_projection=arm_config.feedback_projection,
            critic=arm_config.critic,
            max_concurrency=arm_config.max_concurrency,
        )
        failure_class = ""
        observed_retry_after: float | None = None
        factor_realization_status = "not_observed"
        factor_realization_findings: tuple[str, ...] = ()
        result_sha256 = ""
        result: Any | None = None
        try:
            result = self.runner(
                task=item.case.task,
                provider=ALIBABA_TOKEN_PLAN_PROVIDER,
                provider_config_file=self.provider_config_file,
                secret_file=self.secret_file,
                workspace=item.workspace,
                execution_enabled=False,
                approval_file=None,
                experiment_arm=arm,
                experiment_case=item.case,
                experiment_repeat_index=plan.repeat_index,
                experiment_config=plan.experiment_config,
                approved_molecular_identity=self.source.molecular_identity,
                campaign_preparation_snapshot=(
                    self.campaign_preparation_snapshot
                ),
            )
            _validate_observed_experiment_config(
                result=result,
                expected_config_sha256=plan.experiment_config.config_sha256,
            )
            result_sha256 = _validated_result_digest(result)
            _validate_factor_realization(
                result=result,
                config=plan.experiment_config,
            )
            factor_realization_status = "valid"
            grade = grade_qwen_pyscf_episode(case=item.case, live_result=result)
            failure_class = _classify_live_result(result)
        except Exception as exc:  # noqa: BLE001 - converted to a closed taxonomy
            failure_class, observed_retry_after = _classify_exception(exc)
            if isinstance(exc, _ExperimentFactorRealizationError):
                factor_realization_status = "invalid"
                factor_realization_findings = exc.finding_ids
            failed_result = {
                "terminal_state": "failed",
                "public_transcript": (),
                "successful_tool_calls": 0,
                "failed_tool_calls": 0,
            }
            if not result_sha256:
                result_sha256 = canonical_sha256(
                    {
                        "schema_version": (
                            "chemsmart.controller-failure-result.v1"
                        ),
                        "episode_id": plan.episode_id,
                        "failure_class": failure_class,
                    }
                )
            grade = grade_qwen_pyscf_episode(
                case=item.case, live_result=failed_result
            )
        ledger = _episode_ledger(
            item=item,
            source_sha256=self.source.expected_sha256,
            result_sha256=result_sha256,
            grade=grade,
            failure_class=failure_class,
            factor_realization_status=factor_realization_status,
            factor_realization_findings=factor_realization_findings,
        )
        if self.outcome_observer is not None:
            self.outcome_observer(
                plan=plan,
                case=item.case,
                result=result,
                grade=grade,
                ledger=ledger,
            )
        resume_after = None
        if failure_class in {
            "quota_exhausted",
            "rate_limited",
            "transient_transport",
        }:
            remaining = max(
                0.0,
                (_parse_utc(self.window.deadline_utc) - self._now()).total_seconds(),
            )
            candidate = self.resume_hook(
                episode_id=plan.episode_id,
                failure_class=failure_class,
                observed_retry_after_seconds=observed_retry_after,
                seconds_until_deadline=remaining,
            )
            if candidate is not None:
                candidate = float(candidate)
                if not math.isfinite(candidate) or candidate < 0:
                    raise ContractError("resume hook returned an invalid delay")
                resume_after = candidate
        return _EpisodeOutcome(
            ledger=ledger, resume_after_seconds=resume_after
        )

    def _split_ledger(
        self,
        *,
        split: str,
        started: datetime,
        observed: datetime,
        termination: str,
        ledgers: tuple[QwenPyscfEpisodeLedgerV1, ...],
    ) -> QwenPyscfSplitLedgerV1:
        body = {
            "schema_version": "chemsmart.qwen-split-ledger.v1",
            "campaign_window_sha256": self.window.window_sha256,
            "freeze_manifest_sha256": self.freeze.manifest_sha256,
            "split": split,
            "started_utc": _utc_text(started),
            "observed_at_utc": _utc_text(observed),
            "termination_reason": termination,
            "episode_ledgers": ledgers,
            "last_hypothesis_sha256": (
                ledgers[-1].hypothesis_sha256 if ledgers else ""
            ),
        }
        return QwenPyscfSplitLedgerV1(
            **body, ledger_sha256=canonical_sha256(body)
        )

    def _now(self) -> datetime:
        value = self.clock()
        if value.tzinfo is None or value.utcoffset() is None:
            raise ContractError("campaign clock must return an aware datetime")
        return value.astimezone(timezone.utc)


def _episode_ledger(
    *,
    item: _PreparedEpisode,
    source_sha256: str,
    result_sha256: str,
    grade: QwenPyscfDeterministicGradeV1,
    failure_class: str,
    factor_realization_status: str,
    factor_realization_findings: tuple[str, ...],
) -> QwenPyscfEpisodeLedgerV1:
    safety_violations = set(grade.safety_violations)
    if failure_class == "experiment_integrity_violation":
        safety_violations.add("safety.experiment_config_mismatch")
    body = {
        "schema_version": "chemsmart.qwen-episode-ledger.v1",
        "episode_id": item.plan.episode_id,
        "split": item.case.split,
        "case_sha256": item.case.case_sha256,
        "plan_sha256": item.plan.plan_sha256,
        "hypothesis_sha256": item.plan.hypothesis.hypothesis_sha256,
        "experiment_config_sha256": (
            item.plan.experiment_config.config_sha256
        ),
        "source_sha256": source_sha256,
        "workspace_binding_sha256": item.workspace_binding_sha256,
        "result_sha256": result_sha256,
        "grade_sha256": grade.grade_sha256,
        "session_terminal_state": grade.session_terminal_state,
        "scientific_state": grade.scientific_state,
        "verdict": (
            "inconclusive"
            if factor_realization_status == "invalid"
            else grade.verdict
        ),
        "failure_class": failure_class,
        "factor_realization_status": factor_realization_status,
        "factor_realization_findings": factor_realization_findings,
        "safety_violations": tuple(sorted(safety_violations)),
    }
    return QwenPyscfEpisodeLedgerV1(
        **body, ledger_sha256=canonical_sha256(body)
    )


def _classify_live_result(result: Any) -> str:
    if str(_field(result, "terminal_state", "")) != "failed":
        return ""
    final_text = str(_field(result, "final_text", "")).strip().lower()
    if final_text in _PROVIDER_ERROR_CLASSES:
        return _public_failure_class(final_text)
    if final_text == "provider wall-time budget exhausted":
        return "transient_transport"
    return "episode_failed"


def _classify_exception(exc: Exception) -> tuple[str, float | None]:
    if isinstance(exc, _ExperimentIntegrityError):
        return "experiment_integrity_violation", None
    if isinstance(exc, _ExperimentFactorRealizationError):
        return "experiment_factor_invalid", None
    raw_class = str(getattr(exc, "error_class", "")).strip().lower()
    observed_retry_after = getattr(exc, "retry_after_seconds", None)
    retry_after: float | None = None
    if observed_retry_after is not None:
        try:
            parsed = float(observed_retry_after)
        except (TypeError, ValueError):
            parsed = -1.0
        if math.isfinite(parsed) and parsed >= 0:
            retry_after = parsed
    if raw_class in _PROVIDER_ERROR_CLASSES:
        return _public_failure_class(raw_class), retry_after
    name = type(exc).__name__
    if name in {"DeepSeekProtocolError", "JSONDecodeError"}:
        return "protocol_failure", retry_after
    if isinstance(exc, ContractError):
        return "controller_contract_failure", retry_after
    return "controller_failure", retry_after


def _public_failure_class(raw: str) -> str:
    if raw == "credential_invalid":
        return "credential_revoked"
    if raw == "quota_exhausted":
        return "quota_exhausted"
    if raw == "rate_limited":
        return "rate_limited"
    if raw in _TRANSIENT_ERROR_CLASSES:
        return "transient_transport"
    return "protocol_failure"


def _regular_absolute_file(value: str | Path, label: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute() or not path.is_file() or path.is_symlink():
        raise ContractError(f"{label} must be an absolute regular file")
    return path


def _field(value: Any, name: str, default: Any) -> Any:
    if isinstance(value, Mapping):
        return value.get(name, default)
    return getattr(value, name, default)


def _validated_result_digest(result: Any) -> str:
    observed = require_sha256(
        str(_field(result, "result_sha256", "")), "result_sha256"
    )
    if isinstance(result, Mapping):
        body = {
            key: value for key, value in result.items() if key != "result_sha256"
        }
    elif hasattr(result, "__dict__"):
        body = {
            key: value
            for key, value in result.__dict__.items()
            if key != "result_sha256"
        }
    else:
        raise ContractError("live result has no verifiable public body")
    if observed != canonical_sha256(body):
        raise ContractError("live result digest mismatch")
    return observed


def _validate_observed_experiment_config(
    *, result: Any, expected_config_sha256: str
) -> None:
    expected = require_sha256(
        expected_config_sha256, "expected_config_sha256"
    )
    observations = _field(result, "experiment_observations", {})
    if not isinstance(observations, Mapping):
        raise _ExperimentIntegrityError(
            "live result lacks typed experiment observations"
        )
    if observations.get("experiment_config_sha256") != expected:
        raise _ExperimentIntegrityError(
            "live result used another experiment configuration"
        )
    if observations.get("observed_experiment_config_sha256") != expected:
        raise _ExperimentIntegrityError(
            "live result did not confirm its observed configuration"
        )
    try:
        require_sha256(
            str(observations.get("preparation_sha256") or ""),
            "preparation_sha256",
        )
    except ContractError as exc:
        raise _ExperimentIntegrityError(
            "live result lacks a valid preparation identity"
        ) from exc


def _validate_factor_realization(
    *, result: Any, config: HarnessExperimentConfigV1
) -> None:
    """Require the live session to realize its frozen D, F, and C factors.

    This gate deliberately runs before scientific grading.  A missing worker,
    merge, or critic makes the episode experimentally inadmissible without
    turning a provider/model composition failure into a safety violation.
    """

    observations = _field(result, "experiment_observations", {})
    findings: set[str] = set()
    if not isinstance(observations, Mapping):
        raise _ExperimentFactorRealizationError(
            "experiment.factor.observations_missing"
        )

    gate = observations.get("gate")
    specialists = observations.get("specialists")
    merge = observations.get("merge")
    critic = observations.get("critic")
    provider_sessions = observations.get("provider_sessions")
    specialist_error = str(
        observations.get("nonsecret_specialist_error_class") or ""
    ).strip()
    if observations.get("feedback_projection") != config.feedback_projection:
        findings.add("experiment.factor.feedback_projection_mismatch")

    if not isinstance(gate, Mapping):
        findings.add("experiment.factor.gate_missing")
        gate = {}
    if not isinstance(specialists, (tuple, list)):
        findings.add("experiment.factor.specialists_malformed")
        specialists = ()
    if not isinstance(merge, Mapping):
        findings.add("experiment.factor.merge_missing")
        merge = {}
    if not isinstance(critic, Mapping):
        findings.add("experiment.factor.critic_missing")
        critic = {}
    if not isinstance(provider_sessions, (tuple, list)):
        findings.add("experiment.factor.provider_sessions_malformed")
        provider_sessions = ()

    activated = gate.get("activated") is True
    requested_roles = _factor_roles(
        gate.get("requested_roles"),
        finding_id="experiment.factor.requested_roles_malformed",
        findings=findings,
    )
    accepted_roles = _factor_roles(
        specialists,
        finding_id="experiment.factor.specialist_rows_malformed",
        findings=findings,
        rows=True,
    )
    provider_roles = _factor_roles(
        provider_sessions,
        finding_id="experiment.factor.provider_rows_malformed",
        findings=findings,
        rows=True,
    )
    worker_provider_roles = tuple(
        role for role in provider_roles if role != "critic"
    )
    provider_feedback_modes = tuple(
        str(item.get("feedback_projection") or "")
        for item in provider_sessions
        if isinstance(item, Mapping)
    )
    if any(
        mode != config.feedback_projection for mode in provider_feedback_modes
    ):
        findings.add("experiment.factor.worker_feedback_projection_mismatch")
    for item in provider_sessions:
        if not isinstance(item, Mapping):
            continue
        if not _valid_specialist_output_envelope(
            item.get("output_envelope_receipt")
        ):
            findings.add("experiment.factor.output_envelope_receipt_invalid")

    if config.decomposition and activated:
        if specialist_error:
            findings.add("experiment.factor.specialist_error")
        if accepted_roles != requested_roles:
            findings.add("experiment.factor.specialist_set_incomplete")
        if worker_provider_roles != requested_roles:
            findings.add("experiment.factor.worker_sessions_incomplete")
        if str(merge.get("status") or "") not in {
            "merged",
            "needs_clarification",
        }:
            findings.add("experiment.factor.merge_not_realized")
    else:
        if activated or requested_roles:
            findings.add("experiment.factor.unexpected_gate_activation")
        if accepted_roles or worker_provider_roles or specialist_error:
            findings.add("experiment.factor.unexpected_specialist")
        if str(merge.get("status") or "") != "not_dispatched":
            findings.add("experiment.factor.unexpected_merge")

    critic_status = str(critic.get("status") or "")
    critic_error = str(
        critic.get("nonsecret_error_class") or ""
    ).strip()
    critic_provider_count = provider_roles.count("critic")
    if config.critic:
        if critic_error:
            findings.add("experiment.factor.critic_error")
        if critic_status not in {"complete", "blocked"}:
            findings.add("experiment.factor.critic_not_realized")
        if not _valid_sha256(critic.get("review_sha256")):
            findings.add("experiment.factor.critic_review_missing")
        if not _valid_sha256(critic.get("outcome_sha256")):
            findings.add("experiment.factor.critic_outcome_missing")
        if critic_provider_count != 1:
            findings.add("experiment.factor.critic_session_count")
    else:
        if critic_status != "not_enabled":
            findings.add("experiment.factor.unexpected_critic")
        if critic.get("review_sha256") or critic.get("outcome_sha256"):
            findings.add("experiment.factor.unexpected_critic_receipt")
        if critic_error or critic_provider_count:
            findings.add("experiment.factor.unexpected_critic_session")

    if findings:
        raise _ExperimentFactorRealizationError(*findings)


def _factor_roles(
    value: Any,
    *,
    finding_id: str,
    findings: set[str],
    rows: bool = False,
) -> tuple[str, ...]:
    if not isinstance(value, (tuple, list)):
        findings.add(finding_id)
        return ()
    roles: list[str] = []
    for item in value:
        role = item.get("role") if rows and isinstance(item, Mapping) else item
        if not isinstance(role, str) or not role.strip():
            findings.add(finding_id)
            continue
        roles.append(role.strip())
    if len(roles) != len(set(roles)):
        findings.add(finding_id)
    return tuple(sorted(roles))


def _valid_sha256(value: Any) -> bool:
    try:
        require_sha256(str(value or ""), "factor_receipt_sha256")
    except ContractError:
        return False
    return True


def _valid_specialist_output_envelope(value: Any) -> bool:
    if not isinstance(value, Mapping):
        return False
    expected_keys = {
        "schema_version",
        "mode",
        "raw_text_sha256",
        "normalized_object_sha256",
        "ignored_prefix_bytes",
        "ignored_suffix_bytes",
        "receipt_sha256",
    }
    if set(value) != expected_keys:
        return False
    if value.get("schema_version") != (
        "chemsmart.specialist-output-envelope.v1"
    ):
        return False
    mode = value.get("mode")
    if mode not in {"strict_json", "single_json_object_extracted"}:
        return False
    if not _valid_sha256(value.get("raw_text_sha256")) or not _valid_sha256(
        value.get("normalized_object_sha256")
    ):
        return False
    prefix = value.get("ignored_prefix_bytes")
    suffix = value.get("ignored_suffix_bytes")
    if not isinstance(prefix, int) or not isinstance(suffix, int):
        return False
    if min(prefix, suffix) < 0:
        return False
    if mode == "strict_json" and (prefix or suffix):
        return False
    body = {
        key: item for key, item in value.items() if key != "receipt_sha256"
    }
    return value.get("receipt_sha256") == canonical_sha256(body)


def _experiment_signature(plan: QwenPyscfEpisodePlanV1) -> str:
    """Freeze provider, factors, prompts, tools, order, and budgets.

    ``experiment_id`` is deliberately excluded so a new development repeat
    can carry a unique hypothesis while retaining the preregistered factors.
    """

    config = plan.experiment_config
    return canonical_sha256(
        {
            "schema_version": config.schema_version,
            "decomposition": config.decomposition,
            "feedback_projection": config.feedback_projection,
            "critic": config.critic,
            "provider_id": config.provider_id,
            "model_id": config.model_id,
            "reasoning_effort": config.reasoning_effort,
            "max_concurrency": config.max_concurrency,
            "prompt_sha256": config.prompt_sha256,
            "tool_schema_sha256": config.tool_schema_sha256,
            "task_order_sha256": config.task_order_sha256,
            "token_budget": config.token_budget,
            "tool_call_budget": config.tool_call_budget,
            "wall_time_seconds": config.wall_time_seconds,
        }
    )


def _without_field(value: Any, field_name: str) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.__dict__.items()
        if key != field_name
    }


__all__ = [
    "ApprovedXyzSourceV1",
    "FrozenTransferManifestV1",
    "QwenCampaignWindowV1",
    "QwenPyscfCampaignControllerV1",
    "QwenPyscfEpisodeLedgerV1",
    "QwenPyscfSplitLedgerV1",
    "approved_xyz_source",
    "approved_xyz_source_from_ledger",
    "bind_case_to_approved_xyz",
    "build_frozen_transfer_manifest",
    "build_qwen_campaign_window",
    "prepare_qwen_campaign_plans",
]
