"""Deterministic paired benchmark-plan construction."""

from __future__ import annotations

from typing import Iterable

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.benchmarks.contracts import (
    BenchmarkCaseV1,
    BenchmarkConfigurationV1,
    BenchmarkRunPlanV1,
)


def paired_run_plan(
    cases: Iterable[BenchmarkCaseV1],
    configurations: Iterable[BenchmarkConfigurationV1],
    *,
    repeats: int,
) -> tuple[BenchmarkRunPlanV1, ...]:
    if repeats < 1:
        raise ContractError("paired benchmark requires at least one repeat")
    ordered_cases = tuple(sorted(cases, key=lambda item: item.case_id))
    ordered_configs = tuple(
        sorted(configurations, key=lambda item: item.configuration_id)
    )
    if not ordered_cases or not ordered_configs:
        raise ContractError("paired benchmark requires cases and configurations")
    result = []
    for repeat in range(repeats):
        for case in ordered_cases:
            pairing_key = canonical_sha256(
                {"case_sha256": case.case_sha256, "repeat_index": repeat}
            )
            for config in ordered_configs:
                body = {
                    "run_id": (
                        f"{case.case_id}.{config.configuration_id}.{repeat}"
                    ),
                    "case_sha256": case.case_sha256,
                    "configuration_sha256": config.configuration_sha256,
                    "repeat_index": repeat,
                    "pairing_key": pairing_key,
                }
                result.append(
                    BenchmarkRunPlanV1(
                        **body, run_sha256=canonical_sha256(body)
                    )
                )
    return tuple(result)


__all__ = ["paired_run_plan"]
