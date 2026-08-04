"""Deterministic public-ledger analysis for the Qwen/PySCF campaign.

The analyzer deliberately has no filesystem, provider, or chemistry-engine
access.  It consumes the public split ledgers plus path-free live experiment
observations.  Feedback receipts must be supplied as public event
observations; raw prompts, responses, and private reasoning are neither read
nor needed.

Historical ``tool-feedback-equivalence.v1`` receipts and current
``tool-feedback-projection-receipt.v2`` receipts are separate evidence
families.  Pairing never compares or aggregates their semantic/action
signatures across that version boundary.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import re
from typing import Any, Iterable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    require_sha256,
)
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    QwenPyscfEpisodeLedgerV1,
    QwenPyscfSplitLedgerV1,
)


_EPISODE_ID = re.compile(
    r"^(?P<case>.+)\.d(?P<d>[01])-f(?P<f>full|causal)-c(?P<c>[01])"
    r"\.r(?P<repeat>[0-9]+)$"
)
_OLD_FEEDBACK_SCHEMA = "chemsmart.tool-feedback-equivalence.v1"
_NEW_FEEDBACK_SCHEMA = "chemsmart.tool-feedback-projection-receipt.v2"
_OLD_FAMILY = "historical-semantic-v1"
_NEW_FAMILY = "causal-action-v2"
_H0_ARM_ID = "d0-ffull-c0"
_DFC_ARM_IDS = tuple(
    sorted(
        f"d{decomposition}-f{feedback}-c{critic}"
        for decomposition in (0, 1)
        for feedback in ("full", "causal")
        for critic in (0, 1)
    )
)
_CANDIDATE_RANKING_RULE = (
    "maximize_strict_passes",
    "minimize_strict_failures",
    "minimize_strict_inconclusive",
    "prefer_critic_off_when_scientific_outcomes_tie",
    "minimize_total_tokens",
    "minimize_total_provider_latency_millis",
    "minimize_profile_complexity",
    "canonical_arm_id_tiebreak",
)
_METRIC_NAMES = (
    "model_attempts",
    "input_tokens",
    "output_tokens",
    "reasoning_tokens",
    "provider_latency_millis",
    "successful_tool_calls",
    "failed_tool_calls",
    "tool_calls",
    "canonical_feedback_bytes",
    "provider_feedback_bytes",
)


@dataclass(frozen=True)
class FeedbackBytesByToolV1:
    """Byte accounting for one tool inside one feedback-contract family."""

    receipt_schema_version: str
    contract_family: str
    signature_kind: str
    tool: str
    mode: str
    receipt_count: int
    canonical_bytes: int
    provider_feedback_bytes: int
    bytes_reduced: int

    def __post_init__(self) -> None:
        expected = {
            _OLD_FEEDBACK_SCHEMA: (
                _OLD_FAMILY,
                "semantic_signature_sha256",
            ),
            _NEW_FEEDBACK_SCHEMA: (
                _NEW_FAMILY,
                "causal_action_signature_sha256",
            ),
        }.get(self.receipt_schema_version)
        if expected is None or (
            self.contract_family,
            self.signature_kind,
        ) != expected:
            raise ContractError("feedback byte summary crosses contract families")
        if not self.tool or self.mode not in {"full-v1", "causal-v1"}:
            raise ContractError("feedback byte summary identity is invalid")
        if min(
            self.receipt_count,
            self.canonical_bytes,
            self.provider_feedback_bytes,
        ) < 0 or self.receipt_count == 0:
            raise ContractError("feedback byte summary counters are invalid")
        if self.bytes_reduced != (
            self.canonical_bytes - self.provider_feedback_bytes
        ):
            raise ContractError("feedback byte reduction is inconsistent")


@dataclass(frozen=True)
class QwenPyscfEpisodeAnalysisV1:
    schema_version: str
    episode_id: str
    split: str
    case_id: str
    case_sha256: str
    repeat_index: int
    arm_id: str
    decomposition: bool
    feedback_projection: str
    critic: bool
    ledger_verdict: str
    strict_scientific_verdict: str
    safety_status: str
    safety_violations: tuple[str, ...]
    factor_realization_status: str
    factor_realization_findings: tuple[str, ...]
    failure_class: str
    observation_status: str
    observation_findings: tuple[str, ...]
    decomposition_activated: bool
    model_attempts: int
    input_tokens: int
    output_tokens: int
    reasoning_tokens: int
    provider_latency_millis: int
    successful_tool_calls: int
    failed_tool_calls: int
    specialist_outcomes: tuple[tuple[str, str], ...]
    specialist_merge_status: str
    critic_status: str
    critic_findings: int
    critic_critical_findings: int
    feedback_contract_families: tuple[str, ...]
    feedback_coverage_status: str
    feedback_bytes_by_tool: tuple[FeedbackBytesByToolV1, ...]
    analysis_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-episode-analysis.v1":
            raise ContractError("unsupported Qwen episode analysis")
        require_sha256(self.case_sha256, "case_sha256")
        if self.repeat_index < 0:
            raise ContractError("analysis repeat index must be non-negative")
        if self.ledger_verdict not in {"pass", "fail", "inconclusive"}:
            raise ContractError("analysis ledger verdict is invalid")
        if self.strict_scientific_verdict not in {
            "pass",
            "fail",
            "inconclusive",
        }:
            raise ContractError("strict scientific verdict is invalid")
        if self.safety_status not in {"green", "red"}:
            raise ContractError("analysis safety status is invalid")
        if self.observation_status not in {"valid", "invalid", "missing"}:
            raise ContractError("analysis observation status is invalid")
        if self.feedback_coverage_status not in {
            "complete",
            "incomplete",
            "not_exercised",
        }:
            raise ContractError("feedback coverage status is invalid")
        for values in (
            self.safety_violations,
            self.factor_realization_findings,
            self.observation_findings,
            self.feedback_contract_families,
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError("analysis findings must be canonical")
        if self.specialist_outcomes != tuple(sorted(self.specialist_outcomes)):
            raise ContractError("specialist outcomes must be canonical")
        if min(
            self.model_attempts,
            self.input_tokens,
            self.output_tokens,
            self.reasoning_tokens,
            self.provider_latency_millis,
            self.successful_tool_calls,
            self.failed_tool_calls,
            self.critic_findings,
            self.critic_critical_findings,
        ) < 0:
            raise ContractError("analysis counters must be non-negative")
        body = _without(self, "analysis_sha256")
        if self.analysis_sha256 != canonical_sha256(body):
            raise ContractError("Qwen episode analysis digest mismatch")

    @property
    def tool_calls(self) -> int:
        return self.successful_tool_calls + self.failed_tool_calls

    @property
    def canonical_feedback_bytes(self) -> int:
        return sum(item.canonical_bytes for item in self.feedback_bytes_by_tool)

    @property
    def provider_feedback_bytes(self) -> int:
        return sum(
            item.provider_feedback_bytes for item in self.feedback_bytes_by_tool
        )


@dataclass(frozen=True)
class QwenPyscfPairAnalysisV1:
    schema_version: str
    factor: str
    effect_scope: str
    pair_key: str
    split: str
    case_sha256: str
    repeat_index: int
    held_factors: tuple[tuple[str, str], ...]
    reference_episode_id: str
    candidate_episode_id: str
    status: str
    scientific_comparison: str
    feedback_contract_family: str
    findings: tuple[str, ...]
    metric_deltas: tuple[tuple[str, int], ...]
    provider_feedback_byte_deltas_by_tool: tuple[tuple[str, int], ...]
    analysis_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-pair-analysis.v1":
            raise ContractError("unsupported Qwen pair analysis")
        if self.factor not in {"D", "F", "C"}:
            raise ContractError("paired factor is invalid")
        expected_scope = {
            "D": "decomposition-package-effect",
            "F": "feedback-projection-effect",
            "C": "critic-detection-only",
        }[self.factor]
        if self.effect_scope != expected_scope:
            raise ContractError("paired effect scope is invalid")
        require_sha256(self.case_sha256, "case_sha256")
        if self.status not in {"complete", "inconclusive"}:
            raise ContractError("paired analysis status is invalid")
        if self.scientific_comparison not in {
            "improved",
            "regressed",
            "tied_pass",
            "tied_fail",
            "inconclusive",
        }:
            raise ContractError("paired scientific comparison is invalid")
        if self.findings != tuple(sorted(set(self.findings))):
            raise ContractError("paired findings must be canonical")
        if self.metric_deltas != tuple(sorted(self.metric_deltas)):
            raise ContractError("paired metrics must be canonical")
        if self.provider_feedback_byte_deltas_by_tool != tuple(
            sorted(self.provider_feedback_byte_deltas_by_tool)
        ):
            raise ContractError("paired feedback metrics must be canonical")
        body = _without(self, "analysis_sha256")
        if self.analysis_sha256 != canonical_sha256(body):
            raise ContractError("Qwen pair analysis digest mismatch")


@dataclass(frozen=True)
class QwenPyscfCampaignAnalysisV1:
    schema_version: str
    split_ledger_sha256s: tuple[str, ...]
    episode_analyses: tuple[QwenPyscfEpisodeAnalysisV1, ...]
    pair_analyses: tuple[QwenPyscfPairAnalysisV1, ...]
    strict_passes: int
    strict_failures: int
    strict_inconclusive: int
    safety_red_episodes: int
    factor_invalid_episodes: int
    factor_not_observed_episodes: int
    analysis_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-campaign-analysis.v1":
            raise ContractError("unsupported Qwen campaign analysis")
        if self.split_ledger_sha256s != tuple(
            sorted(set(self.split_ledger_sha256s))
        ):
            raise ContractError("analysis ledgers must be canonical")
        if tuple(item.episode_id for item in self.episode_analyses) != tuple(
            sorted(item.episode_id for item in self.episode_analyses)
        ):
            raise ContractError("episode analyses must be canonical")
        if tuple(item.pair_key for item in self.pair_analyses) != tuple(
            sorted(item.pair_key for item in self.pair_analyses)
        ):
            raise ContractError("pair analyses must be canonical")
        body = _without(self, "analysis_sha256")
        if self.analysis_sha256 != canonical_sha256(body):
            raise ContractError("Qwen campaign analysis digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(asdict(self))


@dataclass(frozen=True)
class QwenPyscfDevelopmentArmSummaryV1:
    """Development-only evidence used to rank one D/F/C arm.

    The preregistered primary matrix has one fresh episode for each distinct
    case in an arm.  Repeated observations are not silently allowed to inflate
    strict pass counts during candidate selection.
    """

    schema_version: str
    arm_id: str
    decomposition: bool
    feedback_projection: str
    critic: bool
    profile_complexity: int
    preregistered_case_count: int
    episode_count: int
    case_ids: tuple[str, ...]
    strict_passes: int
    strict_failures: int
    strict_inconclusive: int
    safety_status: str
    factor_realization_status: str
    observation_status: str
    total_provider_attempts: int
    total_tokens: int
    total_provider_latency_millis: int
    total_tool_failures: int
    eligibility_status: str
    eligibility_findings: tuple[str, ...]
    summary_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != (
            "chemsmart.qwen-development-arm-summary.v1"
        ):
            raise ContractError("unsupported Qwen development arm summary")
        if self.arm_id not in _DFC_ARM_IDS:
            raise ContractError("development arm summary has an unknown arm")
        expected = _arm_factors(self.arm_id)
        if (
            self.decomposition,
            self.feedback_projection,
            self.critic,
        ) != expected:
            raise ContractError("development arm summary factors disagree")
        if self.profile_complexity != sum(
            (
                int(self.decomposition),
                int(self.feedback_projection == "causal-v1"),
                int(self.critic),
            )
        ):
            raise ContractError("development profile complexity is invalid")
        if self.preregistered_case_count < 1 or self.episode_count < 0:
            raise ContractError("development arm counts are invalid")
        if self.case_ids != tuple(sorted(set(self.case_ids))):
            raise ContractError("development arm cases must be canonical")
        if min(
            self.strict_passes,
            self.strict_failures,
            self.strict_inconclusive,
            self.total_provider_attempts,
            self.total_tokens,
            self.total_provider_latency_millis,
            self.total_tool_failures,
        ) < 0:
            raise ContractError("development arm counters must be non-negative")
        if (
            self.strict_passes
            + self.strict_failures
            + self.strict_inconclusive
            != self.episode_count
        ):
            raise ContractError("development strict counts are inconsistent")
        if self.safety_status not in {"green", "red"}:
            raise ContractError("development arm safety status is invalid")
        if self.factor_realization_status not in {"valid", "invalid"}:
            raise ContractError("development factor status is invalid")
        if self.observation_status not in {"valid", "invalid"}:
            raise ContractError("development observation status is invalid")
        if self.eligibility_status not in {"eligible", "ineligible"}:
            raise ContractError("development eligibility status is invalid")
        if self.eligibility_findings != tuple(
            sorted(set(self.eligibility_findings))
        ):
            raise ContractError("development eligibility findings are invalid")
        if (self.eligibility_status == "eligible") != (
            not self.eligibility_findings
        ):
            raise ContractError("development eligibility is inconsistent")
        body = _without(self, "summary_sha256")
        if self.summary_sha256 != canonical_sha256(body):
            raise ContractError("development arm summary digest mismatch")


@dataclass(frozen=True)
class QwenPyscfDevelopmentCandidateSelectionV1:
    """Pure, development-only decision about opening transfer evaluation."""

    schema_version: str
    source_analysis_sha256: str
    preregistered_case_count_per_arm: int
    baseline_arm_id: str
    arm_summaries: tuple[QwenPyscfDevelopmentArmSummaryV1, ...]
    ranking_rule: tuple[str, ...]
    selected_arm_id: str
    selection_status: str
    transfer_status: str
    rationale: tuple[str, ...]
    selection_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != (
            "chemsmart.qwen-development-candidate-selection.v1"
        ):
            raise ContractError("unsupported Qwen candidate selection")
        require_sha256(self.source_analysis_sha256, "source_analysis_sha256")
        if self.preregistered_case_count_per_arm < 1:
            raise ContractError("candidate selection needs a case count")
        if self.baseline_arm_id != _H0_ARM_ID:
            raise ContractError("candidate selection baseline is not H0")
        if tuple(item.arm_id for item in self.arm_summaries) != _DFC_ARM_IDS:
            raise ContractError("candidate selection must summarize all arms")
        if self.ranking_rule != _CANDIDATE_RANKING_RULE:
            raise ContractError("candidate selection ranking rule changed")
        if self.selection_status not in {"selected", "blocked"}:
            raise ContractError("candidate selection status is invalid")
        if self.transfer_status not in {"admissible", "blocked"}:
            raise ContractError("candidate transfer status is invalid")
        selected = {
            item.arm_id: item for item in self.arm_summaries
        }.get(self.selected_arm_id)
        if self.selection_status == "selected":
            if (
                selected is None
                or selected.arm_id == _H0_ARM_ID
                or selected.eligibility_status != "eligible"
                or self.transfer_status != "admissible"
            ):
                raise ContractError("selected development candidate is invalid")
        elif self.selected_arm_id or self.transfer_status != "blocked":
            raise ContractError("blocked selection cannot name a candidate")
        if not self.rationale:
            raise ContractError("candidate selection requires a rationale")
        body = _without(self, "selection_sha256")
        if self.selection_sha256 != canonical_sha256(body):
            raise ContractError("candidate selection digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return canonical_data(asdict(self))


def select_qwen_pyscf_development_candidate(
    *,
    analysis: QwenPyscfCampaignAnalysisV1,
    preregistered_case_count_per_arm: int,
) -> QwenPyscfDevelopmentCandidateSelectionV1:
    """Select a transfer candidate without reading transfer or local state.

    The function accepts a typed in-memory analysis only.  A mixed or
    non-development analysis is rejected before any outcome is ranked so a
    transfer result cannot leak into development selection.
    """

    if not isinstance(analysis, QwenPyscfCampaignAnalysisV1):
        raise ContractError("candidate selection requires a typed analysis")
    if (
        isinstance(preregistered_case_count_per_arm, bool)
        or not isinstance(preregistered_case_count_per_arm, int)
        or preregistered_case_count_per_arm < 1
    ):
        raise ContractError("preregistered case count must be positive")
    if any(item.split != "development" for item in analysis.episode_analyses):
        raise ContractError("candidate selection rejects non-development episodes")
    if any(item.split != "development" for item in analysis.pair_analyses):
        raise ContractError("candidate selection rejects non-development pairs")

    grouped: dict[str, list[QwenPyscfEpisodeAnalysisV1]] = {
        arm_id: [] for arm_id in _DFC_ARM_IDS
    }
    for item in analysis.episode_analyses:
        if item.arm_id not in grouped:
            raise ContractError("candidate selection encountered an unknown arm")
        grouped[item.arm_id].append(item)

    baseline_case_ids = tuple(
        sorted({item.case_id for item in grouped[_H0_ARM_ID]})
    )
    summaries = tuple(
        _summarize_development_arm(
            arm_id=arm_id,
            episodes=tuple(grouped[arm_id]),
            preregistered_case_count=preregistered_case_count_per_arm,
            baseline_case_ids=baseline_case_ids,
        )
        for arm_id in _DFC_ARM_IDS
    )
    by_arm = {item.arm_id: item for item in summaries}
    baseline = by_arm[_H0_ARM_ID]
    eligible = tuple(
        item
        for item in summaries
        if item.arm_id != _H0_ARM_ID
        and item.eligibility_status == "eligible"
    )

    selected_arm_id = ""
    selection_status = "blocked"
    transfer_status = "blocked"
    if baseline.eligibility_status != "eligible":
        rationale = (
            "transfer_blocked: H0 lacks complete safety-green development evidence",
        )
    elif not eligible:
        rationale = (
            "transfer_blocked: no safety-green, factor-valid, "
            "observation-valid non-H0 arm has the preregistered paired cases",
        )
    else:
        selected = min(eligible, key=_candidate_rank_key)
        selected_arm_id = selected.arm_id
        selection_status = "selected"
        transfer_status = "admissible"
        rationale = (
            f"selected {selected.arm_id} by the preregistered lexicographic rule",
            "strict outcomes: "
            f"pass={selected.strict_passes}, fail={selected.strict_failures}, "
            f"inconclusive={selected.strict_inconclusive}",
            "critic-off preference applies only after strict outcomes tie "
            "because C is detection-only",
        )

    body = {
        "schema_version": (
            "chemsmart.qwen-development-candidate-selection.v1"
        ),
        "source_analysis_sha256": analysis.analysis_sha256,
        "preregistered_case_count_per_arm": (
            preregistered_case_count_per_arm
        ),
        "baseline_arm_id": _H0_ARM_ID,
        "arm_summaries": summaries,
        "ranking_rule": _CANDIDATE_RANKING_RULE,
        "selected_arm_id": selected_arm_id,
        "selection_status": selection_status,
        "transfer_status": transfer_status,
        "rationale": rationale,
    }
    return QwenPyscfDevelopmentCandidateSelectionV1(
        **body, selection_sha256=canonical_sha256(body)
    )


def _summarize_development_arm(
    *,
    arm_id: str,
    episodes: tuple[QwenPyscfEpisodeAnalysisV1, ...],
    preregistered_case_count: int,
    baseline_case_ids: tuple[str, ...],
) -> QwenPyscfDevelopmentArmSummaryV1:
    decomposition, feedback_projection, critic = _arm_factors(arm_id)
    case_ids = tuple(sorted({item.case_id for item in episodes}))
    findings: set[str] = set()
    if len(episodes) != preregistered_case_count:
        findings.add("selection.arm.episode_count_mismatch")
    if len(case_ids) != preregistered_case_count:
        findings.add("selection.arm.case_count_mismatch")
    if len(episodes) != len(case_ids):
        findings.add("selection.arm.repeated_case")
    if arm_id != _H0_ARM_ID and case_ids != baseline_case_ids:
        findings.add("selection.arm.h0_case_set_mismatch")
    if any(item.safety_status != "green" for item in episodes):
        findings.add("selection.arm.safety_red")
    if any(item.factor_realization_status != "valid" for item in episodes):
        findings.add("selection.arm.factor_invalid")
    if any(item.observation_status != "valid" for item in episodes):
        findings.add("selection.arm.observation_invalid")

    body = {
        "schema_version": "chemsmart.qwen-development-arm-summary.v1",
        "arm_id": arm_id,
        "decomposition": decomposition,
        "feedback_projection": feedback_projection,
        "critic": critic,
        "profile_complexity": sum(
            (
                int(decomposition),
                int(feedback_projection == "causal-v1"),
                int(critic),
            )
        ),
        "preregistered_case_count": preregistered_case_count,
        "episode_count": len(episodes),
        "case_ids": case_ids,
        "strict_passes": sum(
            item.strict_scientific_verdict == "pass" for item in episodes
        ),
        "strict_failures": sum(
            item.strict_scientific_verdict == "fail" for item in episodes
        ),
        "strict_inconclusive": sum(
            item.strict_scientific_verdict == "inconclusive"
            for item in episodes
        ),
        "safety_status": (
            "green"
            if all(item.safety_status == "green" for item in episodes)
            else "red"
        ),
        "factor_realization_status": (
            "valid"
            if all(
                item.factor_realization_status == "valid"
                for item in episodes
            )
            else "invalid"
        ),
        "observation_status": (
            "valid"
            if all(item.observation_status == "valid" for item in episodes)
            else "invalid"
        ),
        "total_provider_attempts": sum(
            item.model_attempts for item in episodes
        ),
        "total_tokens": sum(
            item.input_tokens + item.output_tokens + item.reasoning_tokens
            for item in episodes
        ),
        "total_provider_latency_millis": sum(
            item.provider_latency_millis for item in episodes
        ),
        "total_tool_failures": sum(
            item.failed_tool_calls for item in episodes
        ),
        "eligibility_status": "eligible" if not findings else "ineligible",
        "eligibility_findings": tuple(sorted(findings)),
    }
    return QwenPyscfDevelopmentArmSummaryV1(
        **body, summary_sha256=canonical_sha256(body)
    )


def _candidate_rank_key(
    item: QwenPyscfDevelopmentArmSummaryV1,
) -> tuple[int | str, ...]:
    return (
        -item.strict_passes,
        item.strict_failures,
        item.strict_inconclusive,
        int(item.critic),
        item.total_tokens,
        item.total_provider_latency_millis,
        item.profile_complexity,
        item.arm_id,
    )


def _arm_factors(arm_id: str) -> tuple[bool, str, bool]:
    match = re.fullmatch(r"d([01])-f(full|causal)-c([01])", arm_id)
    if match is None:
        raise ContractError("candidate selection arm ID is invalid")
    return (
        match.group(1) == "1",
        match.group(2) + "-v1",
        match.group(3) == "1",
    )


def analyze_qwen_pyscf_campaign(
    *,
    split_ledgers: Iterable[QwenPyscfSplitLedgerV1],
    live_observations: Mapping[str, Mapping[str, Any]],
) -> QwenPyscfCampaignAnalysisV1:
    """Analyze public evidence without opening a run directory.

    Each live observation is an object with two keys:

    ``experiment_observations``
        The path-free record returned by ``run_live_agent_session``.
    ``feedback_receipts``
        Direct public feedback receipt objects, or public tool-event payloads
        containing ``feedback_equivalence_receipt``.  The legacy event key is
        accepted because the nested receipt schema determines its meaning.
    """

    ledgers = tuple(split_ledgers)
    episodes: list[QwenPyscfEpisodeLedgerV1] = []
    seen_ids: set[str] = set()
    for split_ledger in ledgers:
        if not isinstance(split_ledger, QwenPyscfSplitLedgerV1):
            raise ContractError("analysis requires typed public split ledgers")
        for episode in split_ledger.episode_ledgers:
            if episode.episode_id in seen_ids:
                raise ContractError("analysis split ledgers repeat an episode")
            seen_ids.add(episode.episode_id)
            episodes.append(episode)
    extra = set(live_observations).difference(seen_ids)
    if extra:
        raise ContractError("live observation has no public episode ledger")

    analyses = tuple(
        sorted(
            (
                _analyze_episode(
                    ledger=episode,
                    observation=live_observations.get(episode.episode_id),
                )
                for episode in episodes
            ),
            key=lambda item: item.episode_id,
        )
    )
    pairs = _build_pair_analyses(analyses)
    body = {
        "schema_version": "chemsmart.qwen-campaign-analysis.v1",
        "split_ledger_sha256s": tuple(
            sorted(item.ledger_sha256 for item in ledgers)
        ),
        "episode_analyses": analyses,
        "pair_analyses": pairs,
        "strict_passes": sum(
            item.strict_scientific_verdict == "pass" for item in analyses
        ),
        "strict_failures": sum(
            item.strict_scientific_verdict == "fail" for item in analyses
        ),
        "strict_inconclusive": sum(
            item.strict_scientific_verdict == "inconclusive"
            for item in analyses
        ),
        "safety_red_episodes": sum(
            item.safety_status == "red" for item in analyses
        ),
        "factor_invalid_episodes": sum(
            item.factor_realization_status == "invalid" for item in analyses
        ),
        "factor_not_observed_episodes": sum(
            item.factor_realization_status == "not_observed"
            for item in analyses
        ),
    }
    return QwenPyscfCampaignAnalysisV1(
        **body, analysis_sha256=canonical_sha256(body)
    )


def _analyze_episode(
    *,
    ledger: QwenPyscfEpisodeLedgerV1,
    observation: Mapping[str, Any] | None,
) -> QwenPyscfEpisodeAnalysisV1:
    arm = _parse_episode_id(ledger.episode_id)
    safety_status = "red" if ledger.safety_violations else "green"
    strict = _strict_verdict(ledger)
    observed = _observe_live_record(
        ledger=ledger,
        expected_feedback=arm["feedback_projection"],
        observation=observation,
    )
    body = {
        "schema_version": "chemsmart.qwen-episode-analysis.v1",
        "episode_id": ledger.episode_id,
        "split": ledger.split,
        "case_id": arm["case_id"],
        "case_sha256": ledger.case_sha256,
        "repeat_index": arm["repeat_index"],
        "arm_id": arm["arm_id"],
        "decomposition": arm["decomposition"],
        "feedback_projection": arm["feedback_projection"],
        "critic": arm["critic"],
        "ledger_verdict": ledger.verdict,
        "strict_scientific_verdict": strict,
        "safety_status": safety_status,
        "safety_violations": ledger.safety_violations,
        "factor_realization_status": ledger.factor_realization_status,
        "factor_realization_findings": ledger.factor_realization_findings,
        "failure_class": ledger.failure_class,
        **observed,
    }
    return QwenPyscfEpisodeAnalysisV1(
        **body, analysis_sha256=canonical_sha256(body)
    )


def _observe_live_record(
    *,
    ledger: QwenPyscfEpisodeLedgerV1,
    expected_feedback: str,
    observation: Mapping[str, Any] | None,
) -> dict[str, Any]:
    empty = {
        "observation_status": "missing",
        "observation_findings": ("analysis.live_observation_missing",),
        "decomposition_activated": False,
        "model_attempts": 0,
        "input_tokens": 0,
        "output_tokens": 0,
        "reasoning_tokens": 0,
        "provider_latency_millis": 0,
        "successful_tool_calls": 0,
        "failed_tool_calls": 0,
        "specialist_outcomes": (),
        "specialist_merge_status": "not_observed",
        "critic_status": "not_observed",
        "critic_findings": 0,
        "critic_critical_findings": 0,
        "feedback_contract_families": (),
        "feedback_coverage_status": "not_exercised",
        "feedback_bytes_by_tool": (),
    }
    if observation is None:
        return empty
    if not isinstance(observation, Mapping):
        return {
            **empty,
            "observation_status": "invalid",
            "observation_findings": ("analysis.live_observation_malformed",),
        }

    findings: set[str] = set()
    experiment = observation.get("experiment_observations")
    if not isinstance(experiment, Mapping):
        return {
            **empty,
            "observation_status": "invalid",
            "observation_findings": (
                "analysis.experiment_observations_missing",
            ),
        }
    if experiment.get("schema_version") != (
        "chemsmart.live-harness-experiment-observations.v1"
    ):
        findings.add("analysis.experiment_observations_schema")
    record_sha = str(experiment.get("record_sha256") or "")
    try:
        require_sha256(record_sha, "record_sha256")
        if record_sha != canonical_sha256(
            {key: value for key, value in experiment.items() if key != "record_sha256"}
        ):
            findings.add("analysis.experiment_observations_digest")
    except ContractError:
        findings.add("analysis.experiment_observations_digest")
    if experiment.get("experiment_config_sha256") != (
        ledger.experiment_config_sha256
    ):
        findings.add("analysis.experiment_config_mismatch")

    usage = experiment.get("usage")
    parsed_usage = _usage(usage, findings=findings)
    specialists = _specialists(experiment.get("specialists"), findings=findings)
    merge = experiment.get("merge")
    if not isinstance(merge, Mapping):
        findings.add("analysis.specialist_merge_malformed")
        merge_status = "not_observed"
    else:
        merge_status = str(merge.get("status") or "not_observed")
    critic = experiment.get("critic")
    critic_status, critic_count, critical_count = _critic(
        critic, findings=findings
    )
    gate = experiment.get("gate")
    if not isinstance(gate, Mapping):
        findings.add("analysis.complexity_gate_malformed")
        activated = False
    else:
        activated = gate.get("activated") is True

    feedback, feedback_families, coverage, feedback_findings = (
        _feedback_summaries(
            observation.get(
                "feedback_receipts",
                experiment.get("feedback_receipts", ()),
            ),
            expected_mode=expected_feedback,
            expected_tool_calls=(
                parsed_usage["successful_tool_calls"]
                + parsed_usage["failed_tool_calls"]
            ),
        )
    )
    findings.update(feedback_findings)
    return {
        "observation_status": "invalid" if findings else "valid",
        "observation_findings": tuple(sorted(findings)),
        "decomposition_activated": activated,
        "model_attempts": parsed_usage["provider_attempts"],
        "input_tokens": parsed_usage["input_tokens"],
        "output_tokens": parsed_usage["output_tokens"],
        "reasoning_tokens": parsed_usage["reasoning_tokens"],
        "provider_latency_millis": parsed_usage["wall_time_millis"],
        "successful_tool_calls": parsed_usage["successful_tool_calls"],
        "failed_tool_calls": parsed_usage["failed_tool_calls"],
        "specialist_outcomes": specialists,
        "specialist_merge_status": merge_status,
        "critic_status": critic_status,
        "critic_findings": critic_count,
        "critic_critical_findings": critical_count,
        "feedback_contract_families": feedback_families,
        "feedback_coverage_status": coverage,
        "feedback_bytes_by_tool": feedback,
    }


def _usage(value: Any, *, findings: set[str]) -> dict[str, int]:
    names = (
        "provider_attempts",
        "input_tokens",
        "output_tokens",
        "reasoning_tokens",
        "wall_time_millis",
        "successful_tool_calls",
        "failed_tool_calls",
    )
    if not isinstance(value, Mapping):
        findings.add("analysis.usage_missing")
        return {name: 0 for name in names}
    result: dict[str, int] = {}
    for name in names:
        raw = value.get(name)
        if isinstance(raw, bool) or not isinstance(raw, int) or raw < 0:
            findings.add("analysis.usage_invalid")
            result[name] = 0
        else:
            result[name] = raw
    return result


def _specialists(
    value: Any, *, findings: set[str]
) -> tuple[tuple[str, str], ...]:
    if not isinstance(value, (tuple, list)):
        findings.add("analysis.specialists_malformed")
        return ()
    rows: list[tuple[str, str]] = []
    for item in value:
        if not isinstance(item, Mapping):
            findings.add("analysis.specialists_malformed")
            continue
        role = str(item.get("role") or "").strip()
        status = str(item.get("status") or "").strip()
        if not role or not status:
            findings.add("analysis.specialists_malformed")
            continue
        rows.append((role, status))
    if len({role for role, _ in rows}) != len(rows):
        findings.add("analysis.specialists_duplicate_role")
    return tuple(sorted(rows))


def _critic(value: Any, *, findings: set[str]) -> tuple[str, int, int]:
    if not isinstance(value, Mapping):
        findings.add("analysis.critic_malformed")
        return "not_observed", 0, 0
    status = str(value.get("status") or "not_observed")
    raw_findings = value.get("findings", ())
    if not isinstance(raw_findings, (tuple, list)):
        findings.add("analysis.critic_findings_malformed")
        return status, 0, 0
    critical = 0
    for item in raw_findings:
        if not isinstance(item, Mapping):
            findings.add("analysis.critic_findings_malformed")
            continue
        if str(item.get("severity") or "").lower() == "critical":
            critical += 1
    return status, len(raw_findings), critical


def _feedback_summaries(
    values: Any,
    *,
    expected_mode: str,
    expected_tool_calls: int,
) -> tuple[
    tuple[FeedbackBytesByToolV1, ...],
    tuple[str, ...],
    str,
    tuple[str, ...],
]:
    findings: set[str] = set()
    if not isinstance(values, (tuple, list)):
        return (), (), "incomplete", ("analysis.feedback_receipts_malformed",)
    groups: dict[tuple[str, str, str, str, str], list[tuple[int, int]]] = {}
    for item in values:
        try:
            receipt = _feedback_receipt(item)
        except ContractError:
            findings.add("analysis.feedback_receipt_invalid")
            continue
        key = (
            receipt["schema"],
            receipt["family"],
            receipt["signature_kind"],
            receipt["tool"],
            receipt["mode"],
        )
        groups.setdefault(key, []).append(
            (receipt["canonical_bytes"], receipt["provider_feedback_bytes"])
        )
        if receipt["mode"] != expected_mode:
            findings.add("analysis.feedback_mode_mismatch")
    summaries = tuple(
        sorted(
            (
                FeedbackBytesByToolV1(
                    receipt_schema_version=key[0],
                    contract_family=key[1],
                    signature_kind=key[2],
                    tool=key[3],
                    mode=key[4],
                    receipt_count=len(rows),
                    canonical_bytes=sum(row[0] for row in rows),
                    provider_feedback_bytes=sum(row[1] for row in rows),
                    bytes_reduced=sum(row[0] - row[1] for row in rows),
                )
                for key, rows in groups.items()
            ),
            key=lambda item: (
                item.receipt_schema_version,
                item.tool,
                item.mode,
            ),
        )
    )
    families = tuple(sorted({item.contract_family for item in summaries}))
    valid_receipts = sum(item.receipt_count for item in summaries)
    if len(families) > 1:
        findings.add("analysis.feedback_contract_mixed_versions")
    if expected_tool_calls == 0 and valid_receipts == 0:
        coverage = "not_exercised"
    elif expected_tool_calls == valid_receipts and not findings:
        coverage = "complete"
    else:
        coverage = "incomplete"
        if expected_tool_calls != valid_receipts:
            findings.add("analysis.feedback_receipt_count_mismatch")
    return summaries, families, coverage, tuple(sorted(findings))


def _feedback_receipt(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ContractError("feedback receipt must be a public object")
    nested = value.get("feedback_equivalence_receipt")
    receipt = nested if isinstance(nested, Mapping) else value
    schema = str(receipt.get("schema_version") or "")
    if schema not in {_OLD_FEEDBACK_SCHEMA, _NEW_FEEDBACK_SCHEMA}:
        raise ContractError("unsupported feedback receipt schema")
    observed_sha = str(receipt.get("receipt_sha256") or "")
    require_sha256(observed_sha, "receipt_sha256")
    if observed_sha != canonical_sha256(
        {key: item for key, item in receipt.items() if key != "receipt_sha256"}
    ):
        raise ContractError("feedback receipt digest mismatch")
    tool = str(receipt.get("tool") or "").strip().lower()
    mode = str(receipt.get("mode") or "").strip().lower()
    canonical_bytes = _nonnegative_int(receipt.get("canonical_bytes"))
    provider_bytes = _nonnegative_int(receipt.get("provider_feedback_bytes"))
    if not tool or mode not in {"full-v1", "causal-v1"}:
        raise ContractError("feedback receipt identity is invalid")
    if canonical_bytes == 0 or provider_bytes == 0:
        raise ContractError("feedback receipt byte counts must be positive")
    if schema == _OLD_FEEDBACK_SCHEMA:
        if receipt.get("verdict") != "causal_semantics_preserved":
            raise ContractError("historical feedback verdict is invalid")
        require_sha256(
            str(receipt.get("semantic_signature_sha256") or ""),
            "semantic_signature_sha256",
        )
        family = _OLD_FAMILY
        signature_kind = "semantic_signature_sha256"
    else:
        expected_contract = (
            "chemsmart.full-tool-result.v1"
            if mode == "full-v1"
            else "chemsmart.causal-action-projection.v1"
        )
        expected_verdict = (
            "full_result_preserved"
            if mode == "full-v1"
            else "causal_action_signature_preserved"
        )
        if (
            receipt.get("projection_contract") != expected_contract
            or receipt.get("verdict") != expected_verdict
        ):
            raise ContractError("current feedback projection contract is invalid")
        require_sha256(
            str(receipt.get("causal_action_signature_sha256") or ""),
            "causal_action_signature_sha256",
        )
        if _nonnegative_int(receipt.get("bytes_reduced"), allow_negative=True) != (
            canonical_bytes - provider_bytes
        ):
            raise ContractError("current feedback byte reduction is invalid")
        family = _NEW_FAMILY
        signature_kind = "causal_action_signature_sha256"
    return {
        "schema": schema,
        "family": family,
        "signature_kind": signature_kind,
        "tool": tool,
        "mode": mode,
        "canonical_bytes": canonical_bytes,
        "provider_feedback_bytes": provider_bytes,
    }


def _nonnegative_int(value: Any, *, allow_negative: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractError("feedback byte count must be an integer")
    if not allow_negative and value < 0:
        raise ContractError("feedback byte count must be non-negative")
    return value


def _build_pair_analyses(
    episodes: tuple[QwenPyscfEpisodeAnalysisV1, ...],
) -> tuple[QwenPyscfPairAnalysisV1, ...]:
    groups: dict[
        tuple[str, str, int],
        dict[tuple[bool, str, bool], QwenPyscfEpisodeAnalysisV1],
    ] = {}
    for episode in episodes:
        key = (episode.split, episode.case_sha256, episode.repeat_index)
        arm = (
            episode.decomposition,
            episode.feedback_projection,
            episode.critic,
        )
        if arm in groups.setdefault(key, {}):
            raise ContractError("analysis contains duplicate D/F/C arm")
        groups[key][arm] = episode

    rows: list[QwenPyscfPairAnalysisV1] = []
    for (split, case_sha, repeat), arms in groups.items():
        d_values = sorted({item[0] for item in arms})
        f_values = sorted({item[1] for item in arms})
        c_values = sorted({item[2] for item in arms})
        for d in d_values:
            for c in c_values:
                rows.append(
                    _pair(
                        factor="F",
                        split=split,
                        case_sha=case_sha,
                        repeat=repeat,
                        held=(("D", str(int(d))), ("C", str(int(c)))),
                        reference=arms.get((d, "full-v1", c)),
                        candidate=arms.get((d, "causal-v1", c)),
                    )
                )
        for f in f_values:
            for c in c_values:
                rows.append(
                    _pair(
                        factor="D",
                        split=split,
                        case_sha=case_sha,
                        repeat=repeat,
                        held=(("F", f), ("C", str(int(c)))),
                        reference=arms.get((False, f, c)),
                        candidate=arms.get((True, f, c)),
                    )
                )
        for d in d_values:
            for f in f_values:
                rows.append(
                    _pair(
                        factor="C",
                        split=split,
                        case_sha=case_sha,
                        repeat=repeat,
                        held=(("D", str(int(d))), ("F", f)),
                        reference=arms.get((d, f, False)),
                        candidate=arms.get((d, f, True)),
                    )
                )
    return tuple(sorted(rows, key=lambda item: item.pair_key))


def _pair(
    *,
    factor: str,
    split: str,
    case_sha: str,
    repeat: int,
    held: tuple[tuple[str, str], ...],
    reference: QwenPyscfEpisodeAnalysisV1 | None,
    candidate: QwenPyscfEpisodeAnalysisV1 | None,
) -> QwenPyscfPairAnalysisV1:
    findings: set[str] = set()
    if reference is None:
        findings.add("analysis.pair.reference_missing")
    if candidate is None:
        findings.add("analysis.pair.candidate_missing")
    if reference is not None and reference.observation_status != "valid":
        findings.add("analysis.pair.reference_observation_invalid")
    if candidate is not None and candidate.observation_status != "valid":
        findings.add("analysis.pair.candidate_observation_invalid")
    if reference is not None and reference.factor_realization_status != "valid":
        findings.add("analysis.pair.reference_factor_invalid")
    if candidate is not None and candidate.factor_realization_status != "valid":
        findings.add("analysis.pair.candidate_factor_invalid")
    if factor == "D" and candidate is not None and not candidate.decomposition_activated:
        findings.add("analysis.pair.decomposition_not_exercised")
    if factor == "F" and reference is not None and candidate is not None:
        if (
            reference.feedback_coverage_status != "complete"
            or candidate.feedback_coverage_status != "complete"
        ):
            findings.add("analysis.pair.feedback_not_exercised")
    family = ""
    if reference is not None and candidate is not None:
        ref_families = reference.feedback_contract_families
        cand_families = candidate.feedback_contract_families
        if ref_families or cand_families:
            if (
                len(ref_families) != 1
                or len(cand_families) != 1
                or ref_families != cand_families
            ):
                findings.add("analysis.pair.feedback_contract_mismatch")
            else:
                family = ref_families[0]

    complete = not findings
    comparison = "inconclusive"
    if (
        complete
        and factor != "C"
        and reference is not None
        and candidate is not None
    ):
        comparison = _scientific_comparison(reference, candidate)
    metric_deltas: tuple[tuple[str, int], ...] = ()
    feedback_deltas: tuple[tuple[str, int], ...] = ()
    if complete and reference is not None and candidate is not None:
        metric_deltas = tuple(
            sorted(
                (
                    name,
                    _metric(candidate, name) - _metric(reference, name),
                )
                for name in _METRIC_NAMES
            )
        )
        ref_bytes = _feedback_by_tool(reference)
        candidate_bytes = _feedback_by_tool(candidate)
        feedback_deltas = tuple(
            sorted(
                (
                    tool,
                    candidate_bytes.get(tool, 0) - ref_bytes.get(tool, 0),
                )
                for tool in set(ref_bytes) | set(candidate_bytes)
            )
        )
    pair_key = canonical_sha256(
        {
            "factor": factor,
            "split": split,
            "case_sha256": case_sha,
            "repeat_index": repeat,
            "held_factors": held,
        }
    )
    body = {
        "schema_version": "chemsmart.qwen-pair-analysis.v1",
        "factor": factor,
        "effect_scope": {
            "D": "decomposition-package-effect",
            "F": "feedback-projection-effect",
            "C": "critic-detection-only",
        }[factor],
        "pair_key": pair_key,
        "split": split,
        "case_sha256": case_sha,
        "repeat_index": repeat,
        "held_factors": tuple(sorted(held)),
        "reference_episode_id": reference.episode_id if reference else "",
        "candidate_episode_id": candidate.episode_id if candidate else "",
        "status": "complete" if complete else "inconclusive",
        "scientific_comparison": comparison,
        "feedback_contract_family": family,
        "findings": tuple(sorted(findings)),
        "metric_deltas": metric_deltas,
        "provider_feedback_byte_deltas_by_tool": feedback_deltas,
    }
    return QwenPyscfPairAnalysisV1(
        **body, analysis_sha256=canonical_sha256(body)
    )


def _parse_episode_id(value: str) -> dict[str, Any]:
    match = _EPISODE_ID.fullmatch(value)
    if match is None:
        raise ContractError("episode ID does not encode a D/F/C arm")
    decomposition = match.group("d") == "1"
    feedback = match.group("f") + "-v1"
    critic = match.group("c") == "1"
    return {
        "case_id": match.group("case"),
        "repeat_index": int(match.group("repeat")),
        "arm_id": (
            f"d{int(decomposition)}-f{match.group('f')}-c{int(critic)}"
        ),
        "decomposition": decomposition,
        "feedback_projection": feedback,
        "critic": critic,
    }


def _strict_verdict(ledger: QwenPyscfEpisodeLedgerV1) -> str:
    if ledger.safety_violations:
        return "fail"
    if ledger.factor_realization_status != "valid":
        return "inconclusive"
    return ledger.verdict


def _scientific_comparison(
    reference: QwenPyscfEpisodeAnalysisV1,
    candidate: QwenPyscfEpisodeAnalysisV1,
) -> str:
    ref = reference.strict_scientific_verdict
    cand = candidate.strict_scientific_verdict
    if "inconclusive" in {ref, cand}:
        return "inconclusive"
    if ref == cand:
        return "tied_pass" if ref == "pass" else "tied_fail"
    return "improved" if cand == "pass" else "regressed"


def _metric(value: QwenPyscfEpisodeAnalysisV1, name: str) -> int:
    if name == "tool_calls":
        return value.tool_calls
    if name == "canonical_feedback_bytes":
        return value.canonical_feedback_bytes
    if name == "provider_feedback_bytes":
        return value.provider_feedback_bytes
    return int(getattr(value, name))


def _feedback_by_tool(
    value: QwenPyscfEpisodeAnalysisV1,
) -> dict[str, int]:
    result: dict[str, int] = {}
    for item in value.feedback_bytes_by_tool:
        result[item.tool] = result.get(item.tool, 0) + item.provider_feedback_bytes
    return result


def _without(value: Any, field_name: str) -> dict[str, Any]:
    return {
        key: item
        for key, item in value.__dict__.items()
        if key != field_name
    }


__all__ = [
    "FeedbackBytesByToolV1",
    "QwenPyscfCampaignAnalysisV1",
    "QwenPyscfDevelopmentArmSummaryV1",
    "QwenPyscfDevelopmentCandidateSelectionV1",
    "QwenPyscfEpisodeAnalysisV1",
    "QwenPyscfPairAnalysisV1",
    "analyze_qwen_pyscf_campaign",
    "select_qwen_pyscf_development_candidate",
]
