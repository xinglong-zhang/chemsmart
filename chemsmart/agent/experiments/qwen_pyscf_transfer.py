"""Development-gated opening of the frozen Qwen/PySCF transfer split."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import re
from typing import Iterable

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    require_sha256,
)
from chemsmart.agent.experiments.qwen_pyscf_analysis import (
    QwenPyscfCampaignAnalysisV1,
    QwenPyscfDevelopmentCandidateSelectionV1,
    select_qwen_pyscf_development_candidate,
)
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    FrozenTransferManifestV1,
    QwenPyscfEpisodeLedgerV1,
    QwenPyscfSplitLedgerV1,
)
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    QwenPyscfCaseSpecV1,
    QwenPyscfEpisodePlanV1,
)


_EPISODE_ID = re.compile(
    r"^(?P<case>.+)\.d(?P<d>[01])-f(?P<f>full|causal)-c(?P<c>[01])"
    r"\.r(?P<repeat>[0-9]+)$"
)
_ARM_ID = re.compile(r"^d[01]-f(?:full|causal)-c[01]$")
_BASELINE_ARM_ID = "d0-ffull-c0"
_DEVELOPMENT_CASES_PER_ARM = 3


@dataclass(frozen=True)
class QwenPyscfTransferPairEvaluationV1:
    """One frozen H0/candidate comparison in the transfer split."""

    schema_version: str
    case_id: str
    repeat_index: int
    baseline_episode_id: str
    candidate_episode_id: str
    baseline_plan_sha256: str
    candidate_plan_sha256: str
    baseline_ledger_sha256: str
    candidate_ledger_sha256: str
    baseline_strict_verdict: str
    candidate_strict_verdict: str
    comparison: str
    pair_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-transfer-pair-evaluation.v1":
            raise ContractError("unsupported Qwen transfer pair evaluation")
        if not self.case_id or self.repeat_index < 0:
            raise ContractError("transfer pair identity is invalid")
        for name in (
            "baseline_plan_sha256",
            "candidate_plan_sha256",
            "baseline_ledger_sha256",
            "candidate_ledger_sha256",
        ):
            require_sha256(getattr(self, name), name)
        allowed_verdicts = {"pass", "fail", "inconclusive"}
        if (
            self.baseline_strict_verdict not in allowed_verdicts
            or self.candidate_strict_verdict not in allowed_verdicts
        ):
            raise ContractError("transfer pair verdict is invalid")
        if self.comparison not in {
            "candidate_better",
            "baseline_better",
            "tied_pass",
            "tied_fail",
            "inconclusive",
        }:
            raise ContractError("transfer pair comparison is invalid")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "pair_sha256"
        }
        if self.pair_sha256 != canonical_sha256(body):
            raise ContractError("transfer pair evaluation digest mismatch")


@dataclass(frozen=True)
class QwenPyscfTransferEvaluationReceiptV1:
    """Terminal evidence binding for the frozen H0/candidate transfer study.

    ``receipt_status`` describes evidence assembly only.  It is deliberately
    not a scientific-success or profile-adoption decision.
    """

    schema_version: str
    implementation_manifest_sha256: str
    implementation_tree_sha256: str
    transfer_freeze_manifest_sha256: str
    development_analysis_sha256: str
    candidate_selection_sha256: str
    transfer_opening_sha256: str
    transfer_analysis_sha256: str
    baseline_arm_id: str
    candidate_arm_id: str
    selected_plan_sha256s: tuple[str, ...]
    split_ledger_sha256s: tuple[str, ...]
    pair_evaluations: tuple[QwenPyscfTransferPairEvaluationV1, ...]
    baseline_strict_counts: tuple[tuple[str, int], ...]
    candidate_strict_counts: tuple[tuple[str, int], ...]
    baseline_safety_red_episodes: int
    candidate_safety_red_episodes: int
    baseline_factor_invalid_episodes: int
    candidate_factor_invalid_episodes: int
    receipt_status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-transfer-evaluation-receipt.v1":
            raise ContractError("unsupported Qwen transfer evaluation receipt")
        for name in (
            "implementation_manifest_sha256",
            "implementation_tree_sha256",
            "transfer_freeze_manifest_sha256",
            "development_analysis_sha256",
            "candidate_selection_sha256",
            "transfer_opening_sha256",
            "transfer_analysis_sha256",
        ):
            require_sha256(getattr(self, name), name)
        if self.baseline_arm_id != _BASELINE_ARM_ID:
            raise ContractError("transfer evaluation baseline is not H0")
        if (
            _ARM_ID.fullmatch(self.candidate_arm_id) is None
            or self.candidate_arm_id == self.baseline_arm_id
        ):
            raise ContractError("transfer evaluation candidate is invalid")
        if self.selected_plan_sha256s != tuple(
            sorted(set(self.selected_plan_sha256s))
        ) or len(self.selected_plan_sha256s) != 80:
            raise ContractError("transfer evaluation requires 80 frozen plans")
        if self.split_ledger_sha256s != tuple(
            sorted(set(self.split_ledger_sha256s))
        ) or not self.split_ledger_sha256s:
            raise ContractError("transfer evaluation ledgers are not canonical")
        for digest in (*self.selected_plan_sha256s, *self.split_ledger_sha256s):
            require_sha256(digest, "transfer evaluation digest")
        if len(self.pair_evaluations) != 40:
            raise ContractError("transfer evaluation requires 40 pairs")
        pair_ids = tuple(
            (item.case_id, item.repeat_index) for item in self.pair_evaluations
        )
        if pair_ids != tuple(sorted(set(pair_ids))):
            raise ContractError("transfer evaluation pairs are not canonical")
        expected_counts = ("fail", "inconclusive", "pass")
        for counts in (
            self.baseline_strict_counts,
            self.candidate_strict_counts,
        ):
            if tuple(name for name, _ in counts) != expected_counts:
                raise ContractError("transfer strict counts are not canonical")
            if any(value < 0 for _, value in counts) or sum(
                value for _, value in counts
            ) != 40:
                raise ContractError("transfer strict counts are inconsistent")
        if min(
            self.baseline_safety_red_episodes,
            self.candidate_safety_red_episodes,
            self.baseline_factor_invalid_episodes,
            self.candidate_factor_invalid_episodes,
        ) < 0:
            raise ContractError("transfer finding counts must be non-negative")
        if self.receipt_status != "complete":
            raise ContractError("transfer evaluation receipt is not complete")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("transfer evaluation receipt digest mismatch")

    def public_record(self) -> dict:
        return canonical_data(asdict(self))


@dataclass(frozen=True)
class QwenPyscfTransferOpeningV1:
    """Exact development-only authorization to expose frozen transfer tasks."""

    schema_version: str
    implementation_manifest_sha256: str
    implementation_tree_sha256: str
    transfer_freeze_manifest_sha256: str
    development_analysis_sha256: str
    candidate_selection_sha256: str
    baseline_arm_id: str
    candidate_arm_id: str
    case_ids: tuple[str, ...]
    repeats_per_case: int
    selected_plan_sha256s: tuple[str, ...]
    paired_schedule: tuple[tuple[str, str], ...]
    provider_exposure_before_selection: int
    comparison_interpretation: str
    opening_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.qwen-transfer-opening.v1":
            raise ContractError("unsupported Qwen transfer opening")
        for name in (
            "implementation_manifest_sha256",
            "implementation_tree_sha256",
            "transfer_freeze_manifest_sha256",
            "development_analysis_sha256",
            "candidate_selection_sha256",
        ):
            require_sha256(getattr(self, name), name)
        if self.baseline_arm_id != _BASELINE_ARM_ID:
            raise ContractError("transfer opening baseline is not H0")
        if (
            _ARM_ID.fullmatch(self.candidate_arm_id) is None
            or self.candidate_arm_id == self.baseline_arm_id
        ):
            raise ContractError("transfer opening requires one non-H0 candidate")
        if (
            self.case_ids != tuple(sorted(set(self.case_ids)))
            or len(self.case_ids) != 4
        ):
            raise ContractError("transfer opening requires four canonical cases")
        if self.repeats_per_case != 10:
            raise ContractError("transfer opening requires ten repeats per case")
        if self.selected_plan_sha256s != tuple(
            sorted(set(self.selected_plan_sha256s))
        ) or len(self.selected_plan_sha256s) != 80:
            raise ContractError("transfer opening requires 80 unique frozen plans")
        for digest in self.selected_plan_sha256s:
            require_sha256(digest, "selected_plan_sha256")
        if any(len(pair) != 2 for pair in self.paired_schedule):
            raise ContractError("transfer opening requires two-episode pairs")
        flattened = tuple(item for pair in self.paired_schedule for item in pair)
        if len(self.paired_schedule) != 40 or len(set(flattened)) != 80:
            raise ContractError("transfer opening schedule is incomplete")
        observed_pairs: set[tuple[str, int]] = set()
        for pair in self.paired_schedule:
            baseline = _parse_episode_identity(pair[0])
            candidate = _parse_episode_identity(pair[1])
            if (
                baseline["arm_id"] != self.baseline_arm_id
                or candidate["arm_id"] != self.candidate_arm_id
            ):
                raise ContractError(
                    "transfer pair order must be H0 then candidate"
                )
            if (
                baseline["case_id"] != candidate["case_id"]
                or baseline["repeat_index"] != candidate["repeat_index"]
            ):
                raise ContractError(
                    "transfer pair must share one case and repeat"
                )
            identity = (
                str(baseline["case_id"]),
                int(baseline["repeat_index"]),
            )
            if identity[0] not in self.case_ids or identity in observed_pairs:
                raise ContractError("transfer pair identity is not canonical")
            observed_pairs.add(identity)
        expected_pairs = {
            (case_id, repeat_index)
            for case_id in self.case_ids
            for repeat_index in range(self.repeats_per_case)
        }
        if observed_pairs != expected_pairs:
            raise ContractError("transfer opening case/repeat coverage differs")
        if self.provider_exposure_before_selection != 0:
            raise ContractError("transfer tasks were exposed before selection")
        if self.comparison_interpretation != "whole-profile comparison":
            raise ContractError("transfer comparison interpretation changed")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "opening_sha256"
        }
        if self.opening_sha256 != canonical_sha256(body):
            raise ContractError("transfer opening digest mismatch")

    def public_record(self) -> dict:
        return canonical_data(asdict(self))


def build_qwen_pyscf_transfer_opening(
    *,
    implementation_manifest_sha256: str,
    implementation_tree_sha256: str,
    freeze: FrozenTransferManifestV1,
    development_analysis: QwenPyscfCampaignAnalysisV1,
    selection: QwenPyscfDevelopmentCandidateSelectionV1,
    provider_exposure_before_selection: int,
    transfer_cases: Iterable[QwenPyscfCaseSpecV1],
    transfer_plans: Iterable[QwenPyscfEpisodePlanV1],
) -> tuple[
    QwenPyscfTransferOpeningV1,
    tuple[tuple[QwenPyscfEpisodePlanV1, QwenPyscfEpisodePlanV1], ...],
]:
    """Open exactly H0/candidate pairs without consulting transfer outcomes."""

    if selection.source_analysis_sha256 != development_analysis.analysis_sha256:
        raise ContractError("candidate selection belongs to another analysis")
    expected_selection = select_qwen_pyscf_development_candidate(
        analysis=development_analysis,
        preregistered_case_count_per_arm=(
            selection.preregistered_case_count_per_arm
        ),
    )
    if expected_selection.selection_sha256 != selection.selection_sha256:
        raise ContractError("candidate selection is not the deterministic result")
    if (
        selection.selection_status != "selected"
        or selection.transfer_status != "admissible"
    ):
        raise ContractError("development evidence did not admit transfer")
    if (
        isinstance(provider_exposure_before_selection, bool)
        or not isinstance(provider_exposure_before_selection, int)
        or provider_exposure_before_selection != 0
    ):
        raise ContractError("transfer preparation exposed a provider")
    development_case_sha256s = tuple(
        sorted(
            {
                item.case_sha256
                for item in development_analysis.episode_analyses
            }
        )
    )
    if (
        selection.preregistered_case_count_per_arm
        != _DEVELOPMENT_CASES_PER_ARM
        or len(freeze.development_case_sha256s)
        != _DEVELOPMENT_CASES_PER_ARM
    ):
        raise ContractError("transfer opening requires three development cases")
    if development_case_sha256s != freeze.development_case_sha256s:
        raise ContractError(
            "development analysis case registry differs from the freeze"
        )
    cases = tuple(transfer_cases)
    case_ids = tuple(sorted(item.case_id for item in cases))
    if len(cases) != 4 or len(set(case_ids)) != 4:
        raise ContractError("transfer opening requires four unique cases")
    if any(item.split != "transfer" for item in cases):
        raise ContractError("transfer opening received a development case")
    case_by_sha = {item.case_sha256: item for item in cases}
    if tuple(sorted(case_by_sha)) != freeze.transfer_case_sha256s:
        raise ContractError("transfer case registry differs from the freeze")

    candidate_arm_id = selection.selected_arm_id
    if _ARM_ID.fullmatch(candidate_arm_id) is None:
        raise ContractError("selected transfer candidate arm is invalid")
    plan_rows = tuple(transfer_plans)
    if len({item.plan_sha256 for item in plan_rows}) != len(plan_rows):
        raise ContractError("transfer plan input repeats a plan digest")
    if len({item.episode_id for item in plan_rows}) != len(plan_rows):
        raise ContractError("transfer plan input repeats an episode")
    selected = []
    for plan in plan_rows:
        parsed = _parse_plan_identity(plan)
        if plan.plan_sha256 not in freeze.transfer_plan_sha256s:
            raise ContractError("transfer plan input contains an unfrozen plan")
        case = case_by_sha.get(plan.case_sha256)
        if case is None:
            raise ContractError("transfer plan references an unknown case")
        if parsed["case_id"] != case.case_id:
            raise ContractError("transfer plan case identity differs")
        _validate_plan_factor_and_pairing(plan=plan, parsed=parsed)
        if parsed["arm_id"] in {_BASELINE_ARM_ID, candidate_arm_id}:
            selected.append((plan, parsed))

    selected_plans = tuple(item[0] for item in selected)
    if len(selected_plans) != 80:
        raise ContractError("selected transfer profile does not have 80 plans")
    if len({item.order_index for item in selected_plans}) != 80:
        raise ContractError("selected transfer plans repeat an order index")
    for case in cases:
        for arm_id in (_BASELINE_ARM_ID, candidate_arm_id):
            repeats = tuple(
                sorted(
                    int(parsed["repeat_index"])
                    for plan, parsed in selected
                    if plan.case_sha256 == case.case_sha256
                    and parsed["arm_id"] == arm_id
                )
            )
            if repeats != tuple(range(10)):
                raise ContractError("transfer case/arm repeats are incomplete")

    grouped: dict[str, list[QwenPyscfEpisodePlanV1]] = {}
    for plan, _ in selected:
        grouped.setdefault(plan.pairing_key, []).append(plan)
    if len(grouped) != 40:
        raise ContractError("transfer plans do not form 40 frozen pairs")
    batches = []
    for rows in grouped.values():
        by_arm = {
            str(_parse_plan_identity(item)["arm_id"]): item for item in rows
        }
        if len(rows) != 2 or set(by_arm) != {
            _BASELINE_ARM_ID,
            candidate_arm_id,
        }:
            raise ContractError("transfer pair does not contain H0 and candidate")
        ordered = (
            by_arm[_BASELINE_ARM_ID],
            by_arm[candidate_arm_id],
        )
        identities = tuple(_parse_plan_identity(item) for item in ordered)
        if (
            identities[0]["case_id"] != identities[1]["case_id"]
            or identities[0]["repeat_index"]
            != identities[1]["repeat_index"]
        ):
            # This branch is unreachable for normal builder output, but keeps
            # a frozen pairing key from joining different cases or repeats.
            raise ContractError("transfer frozen pair identity differs")
        batches.append(ordered)
    paired = tuple(
        rows
        for rows in sorted(
            batches,
            key=lambda items: min(item.order_index for item in items),
        )
    )
    body = {
        "schema_version": "chemsmart.qwen-transfer-opening.v1",
        "implementation_manifest_sha256": require_sha256(
            implementation_manifest_sha256,
            "implementation_manifest_sha256",
        ),
        "implementation_tree_sha256": require_sha256(
            implementation_tree_sha256,
            "implementation_tree_sha256",
        ),
        "transfer_freeze_manifest_sha256": freeze.manifest_sha256,
        "development_analysis_sha256": development_analysis.analysis_sha256,
        "candidate_selection_sha256": selection.selection_sha256,
        "baseline_arm_id": _BASELINE_ARM_ID,
        "candidate_arm_id": candidate_arm_id,
        "case_ids": case_ids,
        "repeats_per_case": 10,
        "selected_plan_sha256s": tuple(
            sorted(item.plan_sha256 for item in selected_plans)
        ),
        "paired_schedule": tuple(
            tuple(item.episode_id for item in pair) for pair in paired
        ),
        "provider_exposure_before_selection": (
            provider_exposure_before_selection
        ),
        "comparison_interpretation": "whole-profile comparison",
    }
    opening = QwenPyscfTransferOpeningV1(
        **body,
        opening_sha256=canonical_sha256(body),
    )
    return opening, paired


def build_qwen_pyscf_transfer_evaluation_receipt(
    *,
    implementation_manifest_sha256: str,
    implementation_tree_sha256: str,
    freeze: FrozenTransferManifestV1,
    development_analysis: QwenPyscfCampaignAnalysisV1,
    selection: QwenPyscfDevelopmentCandidateSelectionV1,
    opening: QwenPyscfTransferOpeningV1,
    transfer_split_ledgers: Iterable[QwenPyscfSplitLedgerV1],
    transfer_analysis: QwenPyscfCampaignAnalysisV1,
) -> QwenPyscfTransferEvaluationReceiptV1:
    """Bind the complete transfer result to its frozen causal ancestry."""

    manifest_sha256 = require_sha256(
        implementation_manifest_sha256,
        "implementation_manifest_sha256",
    )
    tree_sha256 = require_sha256(
        implementation_tree_sha256,
        "implementation_tree_sha256",
    )
    if (
        opening.implementation_manifest_sha256 != manifest_sha256
        or opening.implementation_tree_sha256 != tree_sha256
    ):
        raise ContractError("transfer evaluation implementation differs")
    if opening.transfer_freeze_manifest_sha256 != freeze.manifest_sha256:
        raise ContractError("transfer evaluation freeze differs")
    if (
        opening.development_analysis_sha256
        != development_analysis.analysis_sha256
        or opening.candidate_selection_sha256 != selection.selection_sha256
        or selection.source_analysis_sha256
        != development_analysis.analysis_sha256
    ):
        raise ContractError("transfer evaluation development ancestry differs")
    if (
        opening.baseline_arm_id != selection.baseline_arm_id
        or opening.candidate_arm_id != selection.selected_arm_id
        or selection.transfer_status != "admissible"
    ):
        raise ContractError("transfer evaluation selected profile differs")

    split_ledgers = tuple(transfer_split_ledgers)
    if not split_ledgers or any(
        not isinstance(item, QwenPyscfSplitLedgerV1)
        or item.split != "transfer"
        or item.campaign_window_sha256 != freeze.campaign_window_sha256
        or item.freeze_manifest_sha256 != freeze.manifest_sha256
        for item in split_ledgers
    ):
        raise ContractError("transfer evaluation requires frozen transfer ledgers")
    split_ledger_sha256s = tuple(
        sorted(item.ledger_sha256 for item in split_ledgers)
    )
    if (
        len(split_ledger_sha256s) != len(set(split_ledger_sha256s))
        or transfer_analysis.split_ledger_sha256s != split_ledger_sha256s
    ):
        raise ContractError("transfer analysis ledger ancestry differs")

    ledgers: dict[str, QwenPyscfEpisodeLedgerV1] = {}
    for split_ledger in split_ledgers:
        for episode in split_ledger.episode_ledgers:
            if episode.episode_id in ledgers:
                raise ContractError("transfer evaluation repeats an episode")
            ledgers[episode.episode_id] = episode
    scheduled_episode_ids = tuple(
        item for pair in opening.paired_schedule for item in pair
    )
    if set(ledgers) != set(scheduled_episode_ids) or len(ledgers) != 80:
        raise ContractError("transfer evaluation episode schedule differs")
    if tuple(sorted(item.plan_sha256 for item in ledgers.values())) != (
        opening.selected_plan_sha256s
    ):
        raise ContractError("transfer evaluation plan set differs")

    analyses = {
        item.episode_id: item for item in transfer_analysis.episode_analyses
    }
    if (
        len(analyses) != 80
        or set(analyses) != set(scheduled_episode_ids)
        or any(item.split != "transfer" for item in analyses.values())
    ):
        raise ContractError("transfer evaluation analysis schedule differs")

    pairs = []
    for baseline_id, candidate_id in sorted(
        opening.paired_schedule,
        key=lambda pair: (
            str(_parse_episode_identity(pair[0])["case_id"]),
            int(_parse_episode_identity(pair[0])["repeat_index"]),
        ),
    ):
        baseline_identity = _parse_episode_identity(baseline_id)
        candidate_identity = _parse_episode_identity(candidate_id)
        baseline_analysis = analyses[baseline_id]
        candidate_analysis = analyses[candidate_id]
        if (
            baseline_identity["arm_id"] != opening.baseline_arm_id
            or candidate_identity["arm_id"] != opening.candidate_arm_id
            or baseline_analysis.arm_id != opening.baseline_arm_id
            or candidate_analysis.arm_id != opening.candidate_arm_id
            or baseline_analysis.case_id != baseline_identity["case_id"]
            or candidate_analysis.case_id != candidate_identity["case_id"]
            or baseline_analysis.repeat_index
            != baseline_identity["repeat_index"]
            or candidate_analysis.repeat_index
            != candidate_identity["repeat_index"]
        ):
            raise ContractError("transfer evaluation pair identity differs")
        baseline_ledger = ledgers[baseline_id]
        candidate_ledger = ledgers[candidate_id]
        pair_body = {
            "schema_version": "chemsmart.qwen-transfer-pair-evaluation.v1",
            "case_id": str(baseline_identity["case_id"]),
            "repeat_index": int(baseline_identity["repeat_index"]),
            "baseline_episode_id": baseline_id,
            "candidate_episode_id": candidate_id,
            "baseline_plan_sha256": baseline_ledger.plan_sha256,
            "candidate_plan_sha256": candidate_ledger.plan_sha256,
            "baseline_ledger_sha256": baseline_ledger.ledger_sha256,
            "candidate_ledger_sha256": candidate_ledger.ledger_sha256,
            "baseline_strict_verdict": (
                baseline_analysis.strict_scientific_verdict
            ),
            "candidate_strict_verdict": (
                candidate_analysis.strict_scientific_verdict
            ),
            "comparison": _whole_profile_comparison(
                baseline_analysis.strict_scientific_verdict,
                candidate_analysis.strict_scientific_verdict,
            ),
        }
        pairs.append(
            QwenPyscfTransferPairEvaluationV1(
                **pair_body,
                pair_sha256=canonical_sha256(pair_body),
            )
        )

    pair_evaluations = tuple(pairs)
    baseline_rows = tuple(
        analyses[item.baseline_episode_id] for item in pair_evaluations
    )
    candidate_rows = tuple(
        analyses[item.candidate_episode_id] for item in pair_evaluations
    )
    body = {
        "schema_version": "chemsmart.qwen-transfer-evaluation-receipt.v1",
        "implementation_manifest_sha256": manifest_sha256,
        "implementation_tree_sha256": tree_sha256,
        "transfer_freeze_manifest_sha256": freeze.manifest_sha256,
        "development_analysis_sha256": development_analysis.analysis_sha256,
        "candidate_selection_sha256": selection.selection_sha256,
        "transfer_opening_sha256": opening.opening_sha256,
        "transfer_analysis_sha256": transfer_analysis.analysis_sha256,
        "baseline_arm_id": opening.baseline_arm_id,
        "candidate_arm_id": opening.candidate_arm_id,
        "selected_plan_sha256s": opening.selected_plan_sha256s,
        "split_ledger_sha256s": split_ledger_sha256s,
        "pair_evaluations": pair_evaluations,
        "baseline_strict_counts": _strict_counts(baseline_rows),
        "candidate_strict_counts": _strict_counts(candidate_rows),
        "baseline_safety_red_episodes": sum(
            item.safety_status == "red" for item in baseline_rows
        ),
        "candidate_safety_red_episodes": sum(
            item.safety_status == "red" for item in candidate_rows
        ),
        "baseline_factor_invalid_episodes": sum(
            item.factor_realization_status != "valid" for item in baseline_rows
        ),
        "candidate_factor_invalid_episodes": sum(
            item.factor_realization_status != "valid" for item in candidate_rows
        ),
        "receipt_status": "complete",
    }
    return QwenPyscfTransferEvaluationReceiptV1(
        **body,
        receipt_sha256=canonical_sha256(body),
    )


def _strict_counts(rows: Iterable) -> tuple[tuple[str, int], ...]:
    values = tuple(rows)
    return tuple(
        (verdict, sum(item.strict_scientific_verdict == verdict for item in values))
        for verdict in ("fail", "inconclusive", "pass")
    )


def _whole_profile_comparison(baseline: str, candidate: str) -> str:
    if baseline == candidate == "pass":
        return "tied_pass"
    if baseline == candidate == "fail":
        return "tied_fail"
    if candidate == "pass":
        return "candidate_better"
    if baseline == "pass":
        return "baseline_better"
    return "inconclusive"


def _parse_plan_identity(plan: QwenPyscfEpisodePlanV1) -> dict[str, str | int]:
    if not isinstance(plan, QwenPyscfEpisodePlanV1):
        raise ContractError("transfer schedule requires typed episode plans")
    parsed = _parse_episode_identity(plan.episode_id)
    repeat_index = int(parsed["repeat_index"])
    if repeat_index != plan.repeat_index:
        raise ContractError("transfer repeat identity differs")
    return parsed


def _parse_episode_identity(episode_id: str) -> dict[str, str | int]:
    match = _EPISODE_ID.fullmatch(str(episode_id))
    if match is None:
        raise ContractError("transfer episode identity is invalid")
    arm_id = (
        f"d{match.group('d')}-f{match.group('f')}-c{match.group('c')}"
    )
    repeat_index = int(match.group("repeat"))
    return {
        "case_id": match.group("case"),
        "arm_id": arm_id,
        "repeat_index": repeat_index,
    }


def _validate_plan_factor_and_pairing(
    *,
    plan: QwenPyscfEpisodePlanV1,
    parsed: dict[str, str | int],
) -> None:
    feedback = (
        "causal"
        if plan.experiment_config.feedback_projection == "causal-v1"
        else "full"
    )
    config_arm_id = (
        f"d{int(plan.experiment_config.decomposition)}-"
        f"f{feedback}-c{int(plan.experiment_config.critic)}"
    )
    arm_id = str(parsed["arm_id"])
    repeat_index = int(parsed["repeat_index"])
    case_id = str(parsed["case_id"])
    if (
        config_arm_id != arm_id
        or plan.experiment_config.experiment_id != plan.episode_id.lower()
        or plan.hypothesis.changed_factor != "D/F/C=" + arm_id
        or plan.hypothesis.comparator_id
        != f"{case_id}.{_BASELINE_ARM_ID}.r{repeat_index}"
    ):
        raise ContractError("transfer plan factor identity differs")
    expected_pairing_key = canonical_sha256(
        {
            "case_sha256": plan.case_sha256,
            "repeat_index": repeat_index,
        }
    )
    if plan.pairing_key != expected_pairing_key:
        raise ContractError("transfer plan pairing key differs")


__all__ = [
    "QwenPyscfTransferEvaluationReceiptV1",
    "QwenPyscfTransferOpeningV1",
    "QwenPyscfTransferPairEvaluationV1",
    "build_qwen_pyscf_transfer_evaluation_receipt",
    "build_qwen_pyscf_transfer_opening",
]
