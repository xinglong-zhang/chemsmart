"""Preregisterable Qwen 3.8 Max D x F x C experiment contracts.

This module owns experiment identity and scheduling only.  It composes the
normal ``agent.yaml`` -> ``UnifiedSessionRunner`` path and never implements a
second provider loop, a chemistry engine launcher, or a model-owned readiness
decision.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, Iterable, Protocol

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_sha256,
)
from chemsmart.agent.adaptive_api_campaign import AdaptiveHypothesisV1
from chemsmart.agent.feedback import FEEDBACK_MODES, FULL_FEEDBACK_V1
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_MODEL,
    ALIBABA_TOKEN_PLAN_PROVIDER,
)
from chemsmart.agent.workflows import HarnessExperimentConfigV1


_SCHEMA = "chemsmart.qwen-pyscf-dfc-campaign.v1"
_SPECIALIST_ROLES = (
    "scientific-specialist",
    "pyscf-specialist",
    "dag-specialist",
)
_CRITIC_ROLE = "fresh-read-only-critic"


@dataclass(frozen=True)
class QwenDfcArmV1:
    """One unbound D/F/C arm; runtime evidence binds the canonical config."""

    schema_version: str
    arm_id: str
    decomposition: bool
    feedback_projection: str
    critic: bool
    max_concurrency: int
    specialist_roles: tuple[str, ...]
    arm_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-dfc-arm.v1":
            raise ContractError("unsupported Qwen D/F/C arm")
        expected_id = (
            f"d{int(self.decomposition)}-"
            f"f{'causal' if self.feedback_projection == 'causal-v1' else 'full'}-"
            f"c{int(self.critic)}"
        )
        if self.arm_id != expected_id:
            raise ContractError("D/F/C arm ID is not canonical")
        if self.feedback_projection not in FEEDBACK_MODES:
            raise ContractError("D/F/C feedback projection is unsupported")
        if not 1 <= self.max_concurrency <= 4:
            raise ContractError("campaign concurrency must be within [1, 4]")
        expected_roles = (
            _SPECIALIST_ROLES if self.decomposition else ()
        ) + ((_CRITIC_ROLE,) if self.critic else ())
        if self.specialist_roles != expected_roles:
            raise ContractError("D/F/C roles do not match enabled factors")
        body = _without_field(self, "arm_sha256")
        if self.arm_sha256 != canonical_sha256(body):
            raise ContractError("D/F/C arm digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return {
            **_without_field(self, "arm_sha256"),
            "arm_sha256": self.arm_sha256,
        }


def build_qwen_dfc_arm(
    *,
    decomposition: bool,
    feedback_projection: str = FULL_FEEDBACK_V1,
    critic: bool,
    max_concurrency: int = 1,
) -> QwenDfcArmV1:
    normalized_feedback = str(feedback_projection).strip().lower()
    body = {
        "schema_version": "chemsmart.qwen-dfc-arm.v1",
        "arm_id": (
            f"d{int(bool(decomposition))}-"
            f"f{'causal' if normalized_feedback == 'causal-v1' else 'full'}-"
            f"c{int(bool(critic))}"
        ),
        "decomposition": bool(decomposition),
        "feedback_projection": normalized_feedback,
        "critic": bool(critic),
        "max_concurrency": int(max_concurrency),
        "specialist_roles": (
            _SPECIALIST_ROLES if decomposition else ()
        ) + ((_CRITIC_ROLE,) if critic else ()),
    }
    return QwenDfcArmV1(
        **body, arm_sha256=canonical_sha256(body)
    )


def all_dfc_arms(
    *, max_concurrency: int = 1
) -> tuple[QwenDfcArmV1, ...]:
    """Return all eight orthogonal arms in canonical order."""

    rows = (
        build_qwen_dfc_arm(
            decomposition=decomposition,
            feedback_projection=feedback,
            critic=critic,
            max_concurrency=max_concurrency,
        )
        for decomposition, feedback, critic in product(
            (False, True), ("full-v1", "causal-v1"), (False, True)
        )
    )
    return tuple(sorted(rows, key=lambda item: item.arm_id))


def bind_harness_experiment_config(
    *,
    arm: QwenDfcArmV1,
    experiment_id: str,
    prompt_sha256: str,
    tool_schema_sha256: str,
    task_order_sha256: str,
    token_budget: int,
    tool_call_budget: int,
    wall_time_seconds: int,
) -> HarnessExperimentConfigV1:
    """Bind one arm to the canonical Runtime V2 experiment contract."""

    body = {
        "schema_version": "chemsmart.harness-experiment-config.v1",
        "experiment_id": str(experiment_id).strip().lower(),
        "decomposition": arm.decomposition,
        "feedback_projection": arm.feedback_projection,
        "critic": arm.critic,
        "provider_id": ALIBABA_TOKEN_PLAN_PROVIDER,
        "model_id": ALIBABA_TOKEN_PLAN_MODEL,
        "reasoning_effort": "xhigh",
        "max_concurrency": arm.max_concurrency,
        "prompt_sha256": prompt_sha256,
        "tool_schema_sha256": tool_schema_sha256,
        "task_order_sha256": task_order_sha256,
        "token_budget": int(token_budget),
        "tool_call_budget": int(tool_call_budget),
        "wall_time_seconds": int(wall_time_seconds),
    }
    return HarnessExperimentConfigV1(
        **body, config_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class QwenExperimentPreparationV1:
    """Provider-free observation of the exact coordinator base contract."""

    schema_version: str
    episode_id: str
    case_sha256: str
    arm_id: str
    arm_sha256: str
    repeat_index: int
    task_spec_sha256: str
    artifact_sha256s: tuple[str, ...]
    provider_profile_sha256: str
    experiment_config: HarnessExperimentConfigV1
    provider_calls: int
    engine_calls: int
    approval_files: int
    preparation_sha256: str
    host_snapshot_sha256: str = ""

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-experiment-preparation.v1":
            raise ContractError("unsupported Qwen experiment preparation")
        if self.repeat_index < 0:
            raise ContractError("preparation repeat index must be non-negative")
        for name in (
            "case_sha256",
            "arm_sha256",
            "task_spec_sha256",
            "provider_profile_sha256",
        ):
            require_sha256(getattr(self, name), name)
        if self.host_snapshot_sha256:
            require_sha256(
                self.host_snapshot_sha256, "host_snapshot_sha256"
            )
        if not self.episode_id or not self.arm_id:
            raise ContractError("preparation identity is required")
        if self.experiment_config.experiment_id != self.episode_id.lower():
            raise ContractError("preparation config uses another episode")
        expected_arm_id = (
            f"d{int(self.experiment_config.decomposition)}-"
            "f"
            + (
                "causal"
                if self.experiment_config.feedback_projection == "causal-v1"
                else "full"
            )
            + "-"
            f"c{int(self.experiment_config.critic)}"
        )
        if self.arm_id != expected_arm_id:
            raise ContractError("preparation config differs from its D/F/C arm")
        if self.artifact_sha256s != tuple(
            sorted(set(self.artifact_sha256s))
        ) or not self.artifact_sha256s:
            raise ContractError("preparation artifacts must be canonical")
        for digest in self.artifact_sha256s:
            require_sha256(digest, "artifact_sha256")
        if self.provider_calls or self.engine_calls or self.approval_files:
            raise ContractError("preparation probe must remain read-only")
        body = _without_field(self, "preparation_sha256")
        if not self.host_snapshot_sha256:
            body.pop("host_snapshot_sha256", None)
        if self.preparation_sha256 != canonical_sha256(body):
            raise ContractError("Qwen preparation digest mismatch")

    def public_record(self) -> dict[str, Any]:
        body = _without_field(self, "preparation_sha256")
        if not self.host_snapshot_sha256:
            body.pop("host_snapshot_sha256", None)
        return {**body, "preparation_sha256": self.preparation_sha256}


def build_qwen_experiment_preparation(
    *,
    case: "QwenPyscfCaseSpecV1",
    arm: QwenDfcArmV1,
    repeat_index: int,
    task_spec_sha256: str,
    artifact_sha256s: Iterable[str],
    provider_profile_sha256: str,
    prompt_sha256: str,
    tool_schema_sha256: str,
    task_order_sha256: str,
    token_budget: int,
    tool_call_budget: int,
    wall_time_seconds: int,
    host_snapshot_sha256: str = "",
) -> QwenExperimentPreparationV1:
    episode_id = f"{case.case_id}.{arm.arm_id}.r{int(repeat_index)}"
    config = bind_harness_experiment_config(
        arm=arm,
        experiment_id=episode_id,
        prompt_sha256=prompt_sha256,
        tool_schema_sha256=tool_schema_sha256,
        task_order_sha256=task_order_sha256,
        token_budget=token_budget,
        tool_call_budget=tool_call_budget,
        wall_time_seconds=wall_time_seconds,
    )
    body = {
        "schema_version": "chemsmart.qwen-experiment-preparation.v1",
        "episode_id": episode_id,
        "case_sha256": case.case_sha256,
        "arm_id": arm.arm_id,
        "arm_sha256": arm.arm_sha256,
        "repeat_index": int(repeat_index),
        "task_spec_sha256": task_spec_sha256,
        "artifact_sha256s": tuple(sorted(set(artifact_sha256s))),
        "provider_profile_sha256": provider_profile_sha256,
        "experiment_config": config,
        "provider_calls": 0,
        "engine_calls": 0,
        "approval_files": 0,
    }
    normalized_snapshot_sha256 = str(host_snapshot_sha256).strip()
    if normalized_snapshot_sha256:
        require_sha256(
            normalized_snapshot_sha256, "host_snapshot_sha256"
        )
        body["host_snapshot_sha256"] = normalized_snapshot_sha256
    return QwenExperimentPreparationV1(
        **body, preparation_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class ComplexityGateInputV1:
    """Typed scientific features that can justify bounded delegation."""

    multi_stage: bool = False
    unresolved_scientific_settings: bool = False
    program_substitution: bool = False
    solvent_semantics: bool = False
    gpu_request: bool = False
    excited_state: bool = False
    conflicting_evidence: bool = False
    producer_artifact_edges: bool = False


@dataclass(frozen=True)
class ComplexityGateReceiptV1:
    schema_version: str
    arm_sha256: str
    activated: bool
    reason_ids: tuple[str, ...]
    requested_roles: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.complexity-gate-receipt.v1":
            raise ContractError("unsupported complexity gate receipt")
        require_sha256(self.arm_sha256, "arm_sha256")
        if self.reason_ids != tuple(sorted(set(self.reason_ids))):
            raise ContractError("complexity reasons are not canonical")
        if self.activated != bool(self.requested_roles):
            raise ContractError("complexity activation lacks bounded roles")
        body = _without_field(self, "receipt_sha256")
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("complexity gate receipt digest mismatch")


def evaluate_complexity_gate(
    arm: QwenDfcArmV1,
    signals: ComplexityGateInputV1,
) -> ComplexityGateReceiptV1:
    reasons = tuple(
        sorted(
            key
            for key, enabled in signals.__dict__.items()
            if bool(enabled)
        )
    )
    activated = bool(arm.decomposition and reasons)
    body = {
        "schema_version": "chemsmart.complexity-gate-receipt.v1",
        "arm_sha256": arm.arm_sha256,
        "activated": activated,
        "reason_ids": reasons if activated else (),
        "requested_roles": _SPECIALIST_ROLES if activated else (),
    }
    return ComplexityGateReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


class SpecialistDispatchInterfaceV1(Protocol):
    """Narrow injection point implemented by the scientific-DAG runtime."""

    def dispatch_specialists(
        self, *, packets: tuple[Any, ...], gate: ComplexityGateReceiptV1
    ) -> tuple[Any, ...]: ...

    def dispatch_critic(self, *, review_packet: Any) -> Any: ...


@dataclass(frozen=True)
class QwenPyscfCaseSpecV1:
    """Answer-free development or transfer task with deterministic oracle ID."""

    schema_version: str
    case_id: str
    split: str
    family: str
    task: str
    expected_observation: str
    deterministic_oracle_id: str
    source_sha256s: tuple[str, ...]
    case_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-pyscf-case.v1":
            raise ContractError("unsupported Qwen/PySCF case schema")
        if self.split not in {"development", "transfer"}:
            raise ContractError("case split must be development or transfer")
        if not all(
            str(value).strip()
            for value in (
                self.case_id,
                self.family,
                self.task,
                self.expected_observation,
                self.deterministic_oracle_id,
            )
        ):
            raise ContractError("case fields must not be empty")
        if self.source_sha256s != tuple(sorted(set(self.source_sha256s))):
            raise ContractError("case sources are not canonical")
        body = _without_field(self, "case_sha256")
        if self.case_sha256 != canonical_sha256(body):
            raise ContractError("Qwen/PySCF case digest mismatch")


def build_case_spec(
    *,
    case_id: str,
    split: str,
    family: str,
    task: str,
    expected_observation: str,
    deterministic_oracle_id: str,
    source_sha256s: Iterable[str] = (),
) -> QwenPyscfCaseSpecV1:
    body = {
        "schema_version": "chemsmart.qwen-pyscf-case.v1",
        "case_id": str(case_id),
        "split": str(split),
        "family": str(family),
        "task": str(task),
        "expected_observation": str(expected_observation),
        "deterministic_oracle_id": str(deterministic_oracle_id),
        "source_sha256s": tuple(sorted(set(source_sha256s))),
    }
    return QwenPyscfCaseSpecV1(
        **body, case_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class QwenPyscfEpisodePlanV1:
    schema_version: str
    episode_id: str
    case_sha256: str
    experiment_config: HarnessExperimentConfigV1
    repeat_index: int
    order_index: int
    pairing_key: str
    hypothesis: AdaptiveHypothesisV1
    engine_calls: int
    plan_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-pyscf-episode-plan.v1":
            raise ContractError("unsupported Qwen/PySCF episode plan")
        if min(self.repeat_index, self.order_index) < 0:
            raise ContractError("episode indexes must be non-negative")
        if self.engine_calls != 0:
            raise ContractError("Qwen/PySCF campaign must remain preview-only")
        if self.hypothesis.hypothesis_id != self.episode_id:
            raise ContractError("episode and hypothesis identity differ")
        if self.hypothesis.configuration_sha256 != (
            self.experiment_config.config_sha256
        ):
            raise ContractError("episode hypothesis uses another configuration")
        if (
            self.experiment_config.provider_id != ALIBABA_TOKEN_PLAN_PROVIDER
            or self.experiment_config.model_id != ALIBABA_TOKEN_PLAN_MODEL
            or self.experiment_config.reasoning_effort != "xhigh"
        ):
            raise ContractError("episode configuration is not Qwen Token Plan")
        body = _without_field(self, "plan_sha256")
        if self.plan_sha256 != canonical_sha256(body):
            raise ContractError("Qwen/PySCF episode digest mismatch")


def build_episode_plans(
    *,
    cases: Iterable[QwenPyscfCaseSpecV1],
    arms: Iterable[QwenDfcArmV1],
    repeats: int,
    prompt_sha256: str,
    tool_schema_sha256: str,
    task_order_sha256: str,
    token_budget: int,
    tool_call_budget: int,
    wall_time_seconds: int,
) -> tuple[QwenPyscfEpisodePlanV1, ...]:
    """Build metadata-only paired plans without a provider-attempt cap.

    Use :func:`build_episode_plans_from_preparations` for live/frozen runs;
    caller-supplied prompt or tool digests are not an observed live contract.
    """

    if repeats < 1:
        raise ContractError("D/F/C campaign requires at least one repetition")
    for name, digest in (
        ("prompt_sha256", prompt_sha256),
        ("tool_schema_sha256", tool_schema_sha256),
        ("task_order_sha256", task_order_sha256),
    ):
        require_sha256(digest, name)
    case_rows = tuple(sorted(cases, key=lambda item: item.case_id))
    arm_rows = tuple(
        sorted(arms, key=lambda item: item.arm_id)
    )
    if not case_rows or not arm_rows:
        raise ContractError("D/F/C campaign requires cases and arms")
    unordered: list[tuple[tuple[int, str], dict[str, Any]]] = []
    for repeat_index in range(repeats):
        for case in case_rows:
            pairing_key = canonical_sha256(
                {"case_sha256": case.case_sha256, "repeat_index": repeat_index}
            )
            for arm in arm_rows:
                episode_id = (
                    f"{case.case_id}.{arm.arm_id}.r{repeat_index}"
                )
                experiment_config = bind_harness_experiment_config(
                    arm=arm,
                    experiment_id=episode_id,
                    prompt_sha256=prompt_sha256,
                    tool_schema_sha256=tool_schema_sha256,
                    task_order_sha256=task_order_sha256,
                    token_budget=token_budget,
                    tool_call_budget=tool_call_budget,
                    wall_time_seconds=wall_time_seconds,
                )
                hypothesis_body = {
                    "schema_version": "chemsmart.adaptive-hypothesis.v1",
                    "hypothesis_id": episode_id,
                    "changed_factor": "D/F/C=" + arm.arm_id,
                    "comparator_id": (
                        f"{case.case_id}.d0-ffull-c0.r{repeat_index}"
                    ),
                    "expected_outcome": case.expected_observation,
                    "deterministic_oracle_id": case.deterministic_oracle_id,
                    "source_sha256s": tuple(
                        sorted({case.case_sha256, *case.source_sha256s})
                    ),
                    "prompt_sha256": prompt_sha256,
                    "tool_schema_sha256": tool_schema_sha256,
                    "configuration_sha256": experiment_config.config_sha256,
                    "distinct_from_prior": (
                        "Unique case, D/F/C arm, and repetition tuple; repeated "
                        "tasks measure paired stochastic variation."
                    ),
                }
                hypothesis = AdaptiveHypothesisV1(
                    **hypothesis_body,
                    hypothesis_sha256=canonical_sha256(hypothesis_body),
                )
                row = {
                    "schema_version": "chemsmart.qwen-pyscf-episode-plan.v1",
                    "episode_id": episode_id,
                    "case_sha256": case.case_sha256,
                    "experiment_config": experiment_config,
                    "repeat_index": repeat_index,
                    "order_index": -1,
                    "pairing_key": pairing_key,
                    "hypothesis": hypothesis,
                    "engine_calls": 0,
                }
                counterbalance_key = canonical_sha256(
                    {
                        "schema": _SCHEMA,
                        "repeat": repeat_index,
                        "case": case.case_sha256,
                        "config": experiment_config.config_sha256,
                    }
                )
                unordered.append(((repeat_index, counterbalance_key), row))
    result = []
    for order_index, (_, row) in enumerate(sorted(unordered, key=lambda item: item[0])):
        row["order_index"] = order_index
        result.append(
            QwenPyscfEpisodePlanV1(
                **row, plan_sha256=canonical_sha256(row)
            )
        )
    return tuple(result)


def build_episode_plans_from_preparations(
    *,
    cases: Iterable[QwenPyscfCaseSpecV1],
    arms: Iterable[QwenDfcArmV1],
    preparations: Iterable[QwenExperimentPreparationV1],
) -> tuple[QwenPyscfEpisodePlanV1, ...]:
    """Build frozen plans only from provider-free observed contracts."""

    case_by_sha = {item.case_sha256: item for item in cases}
    arm_by_sha = {item.arm_sha256: item for item in arms}
    rows = tuple(preparations)
    if not rows or not case_by_sha or not arm_by_sha:
        raise ContractError("observed planning requires cases, arms, and probes")
    if len({item.episode_id for item in rows}) != len(rows):
        raise ContractError("observed preparations repeat an episode")
    unordered: list[tuple[tuple[int, str], dict[str, Any]]] = []
    for preparation in rows:
        case = case_by_sha.get(preparation.case_sha256)
        arm = arm_by_sha.get(preparation.arm_sha256)
        if case is None or arm is None:
            raise ContractError("preparation references an unfrozen case or arm")
        if preparation.arm_id != arm.arm_id:
            raise ContractError("preparation arm identity differs")
        config = preparation.experiment_config
        if (
            config.decomposition != arm.decomposition
            or config.feedback_projection != arm.feedback_projection
            or config.critic != arm.critic
            or config.max_concurrency != arm.max_concurrency
        ):
            raise ContractError("preparation config differs from its arm")
        pairing_key = canonical_sha256(
            {
                "case_sha256": case.case_sha256,
                "repeat_index": preparation.repeat_index,
            }
        )
        hypothesis_body = {
            "schema_version": "chemsmart.adaptive-hypothesis.v1",
            "hypothesis_id": preparation.episode_id,
            "changed_factor": "D/F/C=" + arm.arm_id,
            "comparator_id": (
                f"{case.case_id}.d0-ffull-c0.r{preparation.repeat_index}"
            ),
            "expected_outcome": case.expected_observation,
            "deterministic_oracle_id": case.deterministic_oracle_id,
            "source_sha256s": tuple(
                sorted(
                    {
                        case.case_sha256,
                        preparation.preparation_sha256,
                        *(
                            (preparation.host_snapshot_sha256,)
                            if preparation.host_snapshot_sha256
                            else ()
                        ),
                        *case.source_sha256s,
                        *preparation.artifact_sha256s,
                    }
                )
            ),
            "prompt_sha256": config.prompt_sha256,
            "tool_schema_sha256": config.tool_schema_sha256,
            "configuration_sha256": config.config_sha256,
            "distinct_from_prior": (
                "Unique observed case, D/F/C arm, and repetition tuple; "
                "the exact live preparation contract was frozen before use."
            ),
        }
        hypothesis = AdaptiveHypothesisV1(
            **hypothesis_body,
            hypothesis_sha256=canonical_sha256(hypothesis_body),
        )
        plan_body = {
            "schema_version": "chemsmart.qwen-pyscf-episode-plan.v1",
            "episode_id": preparation.episode_id,
            "case_sha256": case.case_sha256,
            "experiment_config": config,
            "repeat_index": preparation.repeat_index,
            "order_index": -1,
            "pairing_key": pairing_key,
            "hypothesis": hypothesis,
            "engine_calls": 0,
        }
        order_key = canonical_sha256(
            {
                "schema": _SCHEMA,
                "preparation": preparation.preparation_sha256,
                "case": case.case_sha256,
                "arm": arm.arm_sha256,
            }
        )
        unordered.append(
            ((preparation.repeat_index, order_key), plan_body)
        )
    result = []
    for order_index, (_, body) in enumerate(
        sorted(unordered, key=lambda item: item[0])
    ):
        body["order_index"] = order_index
        result.append(
            QwenPyscfEpisodePlanV1(
                **body, plan_sha256=canonical_sha256(body)
            )
        )
    return tuple(result)


@dataclass(frozen=True)
class QwenPyscfEpisodeReceiptV1:
    schema_version: str
    plan_sha256: str
    terminal_state: str
    public_result_sha256: str
    deterministic_oracle_verdict: str
    metric_values: tuple[tuple[str, float], ...]
    failure_rule_ids: tuple[str, ...]
    safety_violations: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-pyscf-episode-receipt.v1":
            raise ContractError("unsupported Qwen/PySCF episode receipt")
        if self.terminal_state not in {
            "complete",
            "planned",
            "previewed",
            "blocked",
            "failed",
            "waiting_for_approval",
        }:
            raise ContractError("invalid campaign terminal state")
        if self.deterministic_oracle_verdict not in {
            "pass",
            "fail",
            "inconclusive",
        }:
            raise ContractError("invalid deterministic oracle verdict")
        require_sha256(self.plan_sha256, "plan_sha256")
        require_sha256(self.public_result_sha256, "public_result_sha256")
        if self.metric_values != tuple(sorted(self.metric_values)):
            raise ContractError("campaign metrics are not canonical")
        for name, values in (
            ("failure_rule_ids", self.failure_rule_ids),
            ("safety_violations", self.safety_violations),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} must be sorted and unique")
        body = _without_field(self, "receipt_sha256")
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("campaign episode receipt digest mismatch")


def build_episode_receipt(
    *,
    plan: QwenPyscfEpisodePlanV1,
    terminal_state: str,
    public_result_sha256: str,
    deterministic_oracle_verdict: str,
    metric_values: Iterable[tuple[str, float]] = (),
    failure_rule_ids: Iterable[str] = (),
    safety_violations: Iterable[str] = (),
) -> QwenPyscfEpisodeReceiptV1:
    body = {
        "schema_version": "chemsmart.qwen-pyscf-episode-receipt.v1",
        "plan_sha256": plan.plan_sha256,
        "terminal_state": str(terminal_state),
        "public_result_sha256": require_sha256(
            public_result_sha256, "public_result_sha256"
        ),
        "deterministic_oracle_verdict": str(
            deterministic_oracle_verdict
        ),
        "metric_values": tuple(
            sorted((str(name), float(value)) for name, value in metric_values)
        ),
        "failure_rule_ids": tuple(sorted(set(failure_rule_ids))),
        "safety_violations": tuple(sorted(set(safety_violations))),
    }
    return QwenPyscfEpisodeReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


@dataclass(frozen=True)
class QwenPyscfCampaignReceiptV1:
    schema_version: str
    campaign_id: str
    activation_utc: str
    deadline_utc: str
    configuration_sha256s: tuple[str, ...]
    case_sha256s: tuple[str, ...]
    episode_receipt_sha256s: tuple[str, ...]
    transport_attempts_observed: int
    successful_calls_observed: int
    max_concurrency_observed: int
    termination_reason: str
    last_valid_hypothesis_sha256: str
    safety_violations: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-pyscf-campaign-receipt.v1":
            raise ContractError("unsupported Qwen/PySCF campaign receipt")
        if min(
            self.transport_attempts_observed,
            self.successful_calls_observed,
        ) < 0 or self.successful_calls_observed > self.transport_attempts_observed:
            raise ContractError("campaign call observations are inconsistent")
        if not 1 <= self.max_concurrency_observed <= 4:
            raise ContractError("campaign concurrency observation is invalid")
        if self.termination_reason not in {
            "deadline_reached",
            "quota_exhausted",
            "credential_revoked",
            "no_valid_hypothesis",
            "safety_red_line",
            "active",
        }:
            raise ContractError("campaign termination reason is invalid")
        for name, values in (
            ("configuration_sha256s", self.configuration_sha256s),
            ("case_sha256s", self.case_sha256s),
            ("episode_receipt_sha256s", self.episode_receipt_sha256s),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} must be sorted and unique")
            for digest in values:
                require_sha256(digest, name)
        if self.last_valid_hypothesis_sha256:
            require_sha256(
                self.last_valid_hypothesis_sha256,
                "last_valid_hypothesis_sha256",
            )
        body = _without_field(self, "receipt_sha256")
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("campaign receipt digest mismatch")


def build_campaign_receipt(
    *,
    campaign_id: str,
    activation_utc: str,
    deadline_utc: str,
    configurations: Iterable[HarnessExperimentConfigV1],
    cases: Iterable[QwenPyscfCaseSpecV1],
    episode_receipts: Iterable[QwenPyscfEpisodeReceiptV1],
    transport_attempts_observed: int,
    successful_calls_observed: int,
    max_concurrency_observed: int,
    termination_reason: str,
    last_valid_hypothesis_sha256: str = "",
    safety_violations: Iterable[str] = (),
) -> QwenPyscfCampaignReceiptV1:
    body = {
        "schema_version": "chemsmart.qwen-pyscf-campaign-receipt.v1",
        "campaign_id": str(campaign_id),
        "activation_utc": str(activation_utc),
        "deadline_utc": str(deadline_utc),
        "configuration_sha256s": tuple(
            sorted({item.config_sha256 for item in configurations})
        ),
        "case_sha256s": tuple(
            sorted({item.case_sha256 for item in cases})
        ),
        "episode_receipt_sha256s": tuple(
            sorted({item.receipt_sha256 for item in episode_receipts})
        ),
        "transport_attempts_observed": int(transport_attempts_observed),
        "successful_calls_observed": int(successful_calls_observed),
        "max_concurrency_observed": int(max_concurrency_observed),
        "termination_reason": str(termination_reason),
        "last_valid_hypothesis_sha256": str(
            last_valid_hypothesis_sha256
        ),
        "safety_violations": tuple(sorted(set(safety_violations))),
    }
    return QwenPyscfCampaignReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _without_field(value: Any, field_name: str) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.__dict__.items()
        if key != field_name
    }


__all__ = [
    "ComplexityGateInputV1",
    "ComplexityGateReceiptV1",
    "QwenDfcArmV1",
    "QwenExperimentPreparationV1",
    "QwenPyscfCampaignReceiptV1",
    "QwenPyscfCaseSpecV1",
    "QwenPyscfEpisodePlanV1",
    "QwenPyscfEpisodeReceiptV1",
    "SpecialistDispatchInterfaceV1",
    "all_dfc_arms",
    "build_case_spec",
    "build_campaign_receipt",
    "build_episode_plans",
    "build_episode_plans_from_preparations",
    "build_episode_receipt",
    "bind_harness_experiment_config",
    "build_qwen_dfc_arm",
    "build_qwen_experiment_preparation",
    "evaluate_complexity_gate",
]
