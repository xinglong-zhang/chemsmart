"""Reusable benchmark manifests without paper answers or static CLI state."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import ContractError, canonical_sha256


@dataclass(frozen=True)
class BenchmarkOracleV1:
    oracle_id: str
    oracle_version: str
    implementation_sha256: str
    deterministic: bool = True

    def __post_init__(self) -> None:
        if not self.deterministic:
            raise ContractError("primary benchmark oracles must be deterministic")


@dataclass(frozen=True)
class BenchmarkCaseV1:
    case_id: str
    source_artifact_sha256s: tuple[str, ...]
    task_spec_sha256: str
    geometry_artifact_sha256: str
    held_out: bool
    oracle: BenchmarkOracleV1
    case_sha256: str

    def __post_init__(self) -> None:
        body = {
            "case_id": self.case_id,
            "source_artifact_sha256s": self.source_artifact_sha256s,
            "task_spec_sha256": self.task_spec_sha256,
            "geometry_artifact_sha256": self.geometry_artifact_sha256,
            "held_out": self.held_out,
            "oracle": self.oracle,
        }
        if self.case_sha256 != canonical_sha256(body):
            raise ContractError("benchmark case digest mismatch")


@dataclass(frozen=True)
class BenchmarkConfigurationV1:
    configuration_id: str
    factor_levels: tuple[tuple[str, bool], ...]
    provider_model: str
    prompt_sha256: str
    tool_schema_sha256: str
    joined_capability_sha256: str
    budget_sha256: str
    configuration_sha256: str

    def __post_init__(self) -> None:
        names = tuple(item[0] for item in self.factor_levels)
        if names != tuple(sorted(set(names))):
            raise ContractError("benchmark factors must be sorted and unique")
        body = {
            "configuration_id": self.configuration_id,
            "factor_levels": self.factor_levels,
            "provider_model": self.provider_model,
            "prompt_sha256": self.prompt_sha256,
            "tool_schema_sha256": self.tool_schema_sha256,
            "joined_capability_sha256": self.joined_capability_sha256,
            "budget_sha256": self.budget_sha256,
        }
        if self.configuration_sha256 != canonical_sha256(body):
            raise ContractError("benchmark configuration digest mismatch")


@dataclass(frozen=True)
class BenchmarkRunPlanV1:
    run_id: str
    case_sha256: str
    configuration_sha256: str
    repeat_index: int
    pairing_key: str
    run_sha256: str

    def __post_init__(self) -> None:
        if self.repeat_index < 0:
            raise ContractError("repeat index must be non-negative")
        body = {
            "run_id": self.run_id,
            "case_sha256": self.case_sha256,
            "configuration_sha256": self.configuration_sha256,
            "repeat_index": self.repeat_index,
            "pairing_key": self.pairing_key,
        }
        if self.run_sha256 != canonical_sha256(body):
            raise ContractError("benchmark run plan digest mismatch")


@dataclass(frozen=True)
class BenchmarkRunReceiptV1:
    run_sha256: str
    event_stream_sha256: str
    terminal_state: str
    metric_values: tuple[tuple[str, float], ...]
    failure_rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.terminal_state not in {
            "complete",
            "planned",
            "failed",
            "blocked",
            "waiting_for_approval",
        }:
            raise ContractError("benchmark run lacks a terminal state")
        names = tuple(item[0] for item in self.metric_values)
        if names != tuple(sorted(set(names))):
            raise ContractError("benchmark metric names must be sorted and unique")
        body = {
            "run_sha256": self.run_sha256,
            "event_stream_sha256": self.event_stream_sha256,
            "terminal_state": self.terminal_state,
            "metric_values": self.metric_values,
            "failure_rule_ids": self.failure_rule_ids,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("benchmark run receipt digest mismatch")


__all__ = [
    "BenchmarkCaseV1",
    "BenchmarkConfigurationV1",
    "BenchmarkOracleV1",
    "BenchmarkRunPlanV1",
    "BenchmarkRunReceiptV1",
]
