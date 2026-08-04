"""Projection-independent host observations for scientific experiment grading.

Provider feedback is an experimental factor.  It is therefore not a valid
source for the deterministic outcome oracle: ``full-v1`` and ``causal-v1`` may
show the provider different views of the same host action.  This module builds
one path-free grading bundle from the complete canonical tool results retained
in Runtime V2 ``tool_succeeded`` and ``tool_failed`` events.

The bundle deliberately retains typed scientific settings, workflow nodes and
edges, receipt identities, statuses, findings, and preview semantics.  It
never carries local paths, argv or shell text, rendered project/native input,
process output, credentials, or provider reasoning.
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
from chemsmart.agent.runtime.events import EventKind, RuntimeEvent


_HOST_TOOL_OBSERVATION_SCHEMA = "chemsmart.host-tool-observation.v1"
_HOST_ORACLE_BUNDLE_SCHEMA = "chemsmart.host-oracle-input-bundle.v1"

_FORBIDDEN_KEYS = frozenset(
    {
        "argv",
        "canonical_command",
        "cli_value",
        "command",
        "cwd",
        "directory",
        "event_store_path",
        "file",
        "filename",
        "files",
        "generated_input",
        "generated_script",
        "input_text",
        "native_input",
        "native_input_text",
        "path",
        "project_yaml",
        "public_transcript",
        "raw_output",
        "raw_response",
        "rendered_command",
        "rendered_yaml",
        "run_directory",
        "script",
        "script_source",
        "script_text",
        "shell",
        "stderr",
        "stdout",
        "workspace",
        "yaml_text",
    }
)
_SAFE_SEMANTIC_PATH_KEYS = frozenset({"command_path", "field_path"})
_SENSITIVE_KEY_FRAGMENTS = (
    "api_key",
    "authorization",
    "cookie",
    "credential",
    "password",
    "private_reasoning",
    "reasoning_content",
    "secret",
)
_ABSOLUTE_WINDOWS_PATH = re.compile(r"^[A-Za-z]:[\\/]")
_DROP = object()


@dataclass(frozen=True)
class HostToolObservationV1:
    """One canonical host action projected for a deterministic oracle."""

    schema_version: str
    ordinal: int
    tool_name: str
    host_status: str
    canonical_result_sha256: str
    oracle_result: Mapping[str, Any]
    oracle_result_sha256: str
    error_class: str
    rule_ids: tuple[str, ...]
    omitted_fields: tuple[str, ...]
    observation_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _HOST_TOOL_OBSERVATION_SCHEMA:
            raise ContractError("unsupported host tool observation schema")
        if self.ordinal < 1:
            raise ContractError("host tool observation ordinal must be positive")
        if not str(self.tool_name).strip():
            raise ContractError("host tool observation requires a tool name")
        if self.host_status not in {"succeeded", "failed"}:
            raise ContractError("host tool observation status is invalid")
        require_sha256(self.canonical_result_sha256, "canonical_result_sha256")
        if not isinstance(self.oracle_result, Mapping):
            raise ContractError("host oracle result must be an object")
        canonical_result = canonical_data(dict(self.oracle_result))
        _validate_oracle_projection(canonical_result)
        if self.oracle_result_sha256 != canonical_sha256(canonical_result):
            raise ContractError("host oracle result digest mismatch")
        if self.rule_ids != tuple(sorted(set(self.rule_ids))):
            raise ContractError("host observation rule IDs are not canonical")
        if self.omitted_fields != tuple(sorted(set(self.omitted_fields))):
            raise ContractError("host observation omissions are not canonical")
        body = {
            key: value
            for key, value in asdict(self).items()
            if key != "observation_sha256"
        }
        if self.observation_sha256 != canonical_sha256(body):
            raise ContractError("host tool observation digest mismatch")

    @classmethod
    def from_record(cls, value: Mapping[str, Any]) -> "HostToolObservationV1":
        if not isinstance(value, Mapping):
            raise ContractError("host tool observation must be an object")
        result = value.get("oracle_result")
        if not isinstance(result, Mapping):
            raise ContractError("host tool observation lacks an oracle result")
        raw_rule_ids = value.get("rule_ids") or ()
        raw_omitted = value.get("omitted_fields") or ()
        if not isinstance(raw_rule_ids, (tuple, list)) or not isinstance(
            raw_omitted, (tuple, list)
        ):
            raise ContractError("host tool observation sequences are malformed")
        return cls(
            schema_version=str(value.get("schema_version") or ""),
            ordinal=_record_int(value.get("ordinal"), "ordinal"),
            tool_name=str(value.get("tool_name") or ""),
            host_status=str(value.get("host_status") or ""),
            canonical_result_sha256=str(
                value.get("canonical_result_sha256") or ""
            ),
            oracle_result=canonical_data(dict(result)),
            oracle_result_sha256=str(value.get("oracle_result_sha256") or ""),
            error_class=str(value.get("error_class") or ""),
            rule_ids=tuple(str(item) for item in raw_rule_ids),
            omitted_fields=tuple(str(item) for item in raw_omitted),
            observation_sha256=str(value.get("observation_sha256") or ""),
        )

    def public_record(self) -> dict[str, Any]:
        return canonical_data(asdict(self))


@dataclass(frozen=True)
class HostOracleInputBundleV1:
    """F-invariant scientific tool evidence bound to one coordinator stream."""

    schema_version: str
    session_id: str
    event_stream_head_sha256: str
    observations: tuple[HostToolObservationV1, ...]
    successful_tool_calls: int
    failed_tool_calls: int
    tool_counts: tuple[tuple[str, int, int], ...]
    tool_actions_sha256: str
    bundle_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _HOST_ORACLE_BUNDLE_SCHEMA:
            raise ContractError("unsupported host oracle bundle schema")
        if not self.session_id.strip():
            raise ContractError("host oracle bundle requires a session ID")
        require_sha256(
            self.event_stream_head_sha256, "event_stream_head_sha256"
        )
        if min(self.successful_tool_calls, self.failed_tool_calls) < 0:
            raise ContractError("host oracle tool counts must be non-negative")
        if tuple(item.ordinal for item in self.observations) != tuple(
            range(1, len(self.observations) + 1)
        ):
            raise ContractError("host tool observations must be contiguous")
        expected_successful = sum(
            item.host_status == "succeeded" for item in self.observations
        )
        expected_failed = sum(
            item.host_status == "failed" for item in self.observations
        )
        if (
            self.successful_tool_calls != expected_successful
            or self.failed_tool_calls != expected_failed
        ):
            raise ContractError("host oracle counts do not match observations")
        expected_tool_counts = _tool_counts(self.observations)
        if self.tool_counts != expected_tool_counts:
            raise ContractError("host oracle per-tool counts are inconsistent")
        expected_actions = canonical_sha256(
            tuple(item.public_record() for item in self.observations)
        )
        if self.tool_actions_sha256 != expected_actions:
            raise ContractError("host oracle tool-action digest mismatch")
        body = {
            key: value
            for key, value in asdict(self).items()
            if key != "bundle_sha256"
        }
        if self.bundle_sha256 != canonical_sha256(body):
            raise ContractError("host oracle bundle digest mismatch")

    @classmethod
    def from_record(cls, value: Mapping[str, Any]) -> "HostOracleInputBundleV1":
        if not isinstance(value, Mapping):
            raise ContractError("host oracle bundle must be an object")
        raw_observations = value.get("observations")
        if not isinstance(raw_observations, (tuple, list)):
            raise ContractError("host oracle bundle lacks observations")
        raw_counts = value.get("tool_counts")
        if not isinstance(raw_counts, (tuple, list)):
            raise ContractError("host oracle bundle lacks per-tool counts")
        tool_counts: list[tuple[str, int, int]] = []
        for row in raw_counts:
            if not isinstance(row, (tuple, list)) or len(row) != 3:
                raise ContractError("host oracle tool count row is malformed")
            tool_counts.append(
                (
                    str(row[0]),
                    _record_int(row[1], "successful tool count"),
                    _record_int(row[2], "failed tool count"),
                )
            )
        return cls(
            schema_version=str(value.get("schema_version") or ""),
            session_id=str(value.get("session_id") or ""),
            event_stream_head_sha256=str(
                value.get("event_stream_head_sha256") or ""
            ),
            observations=tuple(
                HostToolObservationV1.from_record(item)
                for item in raw_observations
            ),
            successful_tool_calls=_record_int(
                value.get("successful_tool_calls"), "successful_tool_calls"
            ),
            failed_tool_calls=_record_int(
                value.get("failed_tool_calls"), "failed_tool_calls"
            ),
            tool_counts=tuple(tool_counts),
            tool_actions_sha256=str(value.get("tool_actions_sha256") or ""),
            bundle_sha256=str(value.get("bundle_sha256") or ""),
        )

    def public_record(self) -> dict[str, Any]:
        return canonical_data(asdict(self))


def build_host_oracle_input_bundle(
    *,
    events: Iterable[RuntimeEvent],
    session_id: str,
    event_stream_head_sha256: str,
    successful_tool_calls: int,
    failed_tool_calls: int,
) -> HostOracleInputBundleV1:
    """Reconstruct the grading input from full coordinator Runtime V2 events."""

    rows = tuple(events)
    if not rows:
        raise ContractError("host oracle reconstruction requires Runtime V2 events")
    if any(item.session_id != session_id for item in rows):
        raise ContractError("host oracle events contain another session")
    previous_hash = ""
    for expected_sequence, event in enumerate(rows, start=1):
        if event.sequence != expected_sequence or event.previous_hash != previous_hash:
            raise ContractError("host oracle event chain is not contiguous")
        if not event.verify_hash():
            raise ContractError("host oracle event hash is invalid")
        previous_hash = event.event_hash
    if rows[-1].event_hash != event_stream_head_sha256:
        raise ContractError("host oracle stream head differs from the live result")
    observations: list[HostToolObservationV1] = []
    for event in rows:
        if event.kind not in {
            EventKind.TOOL_SUCCEEDED.value,
            EventKind.TOOL_FAILED.value,
        }:
            continue
        observations.append(
            _observation_from_event(event, ordinal=len(observations) + 1)
        )
    body = {
        "schema_version": _HOST_ORACLE_BUNDLE_SCHEMA,
        "session_id": str(session_id),
        "event_stream_head_sha256": require_sha256(
            event_stream_head_sha256, "event_stream_head_sha256"
        ),
        "observations": tuple(observations),
        "successful_tool_calls": int(successful_tool_calls),
        "failed_tool_calls": int(failed_tool_calls),
        "tool_counts": _tool_counts(tuple(observations)),
        "tool_actions_sha256": canonical_sha256(
            tuple(item.public_record() for item in observations)
        ),
    }
    return HostOracleInputBundleV1(
        **body, bundle_sha256=canonical_sha256(body)
    )


def _observation_from_event(
    event: RuntimeEvent, *, ordinal: int
) -> HostToolObservationV1:
    payload = event.payload
    canonical_result = payload.get("canonical_result")
    if not isinstance(canonical_result, Mapping):
        raise ContractError("tool event lacks its complete canonical result")
    canonical_result = canonical_data(dict(canonical_result))
    tool_name = str(payload.get("tool") or "").strip()
    if not tool_name or str(canonical_result.get("tool") or "") != tool_name:
        raise ContractError("tool event and canonical result identity differ")
    canonical_result_sha256 = canonical_sha256(canonical_result)
    observed_result_sha256 = str(payload.get("result_sha256") or "")
    if observed_result_sha256 and observed_result_sha256 != canonical_result_sha256:
        raise ContractError("tool event canonical result digest mismatch")
    nested = canonical_result.get("result")
    if isinstance(nested, Mapping):
        semantic_source: Mapping[str, Any] = nested
    else:
        semantic_source = {
            key: value
            for key, value in canonical_result.items()
            if key not in {"schema_version", "tool", "result"}
        }
    omitted: list[str] = []
    projected = _project_oracle_value(
        semantic_source, path="$.result", omitted=omitted
    )
    if projected is _DROP or not isinstance(projected, Mapping):
        projected = {}
    projected = canonical_data(dict(projected))
    body = {
        "schema_version": _HOST_TOOL_OBSERVATION_SCHEMA,
        "ordinal": int(ordinal),
        "tool_name": tool_name,
        "host_status": (
            "succeeded"
            if event.kind == EventKind.TOOL_SUCCEEDED.value
            else "failed"
        ),
        "canonical_result_sha256": canonical_result_sha256,
        "oracle_result": projected,
        "oracle_result_sha256": canonical_sha256(projected),
        "error_class": str(payload.get("error_class") or ""),
        "rule_ids": tuple(
            sorted({str(item) for item in payload.get("rule_ids") or ()})
        ),
        "omitted_fields": tuple(sorted(set(omitted))),
    }
    return HostToolObservationV1(
        **body, observation_sha256=canonical_sha256(body)
    )


def _tool_counts(
    observations: tuple[HostToolObservationV1, ...],
) -> tuple[tuple[str, int, int], ...]:
    counts: dict[str, list[int]] = {}
    for item in observations:
        row = counts.setdefault(item.tool_name, [0, 0])
        row[0 if item.host_status == "succeeded" else 1] += 1
    return tuple(
        (tool, values[0], values[1])
        for tool, values in sorted(counts.items())
    )


def _record_int(value: Any, field_name: str) -> int:
    try:
        return int(value if value is not None else 0)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"{field_name} must be an integer") from exc


def _project_oracle_value(
    value: Any, *, path: str, omitted: list[str]
) -> Any:
    if isinstance(value, Mapping):
        result: dict[str, Any] = {}
        for key, child in sorted(value.items(), key=lambda pair: str(pair[0])):
            name = str(key)
            child_path = f"{path}.{name}"
            if _forbidden_key(name):
                omitted.append(child_path)
                continue
            projected = _project_oracle_value(
                child, path=child_path, omitted=omitted
            )
            if projected is _DROP:
                omitted.append(child_path)
                continue
            result[name] = projected
        return result
    if isinstance(value, (tuple, list)):
        result = []
        for index, child in enumerate(value):
            child_path = f"{path}[{index}]"
            projected = _project_oracle_value(
                child, path=child_path, omitted=omitted
            )
            if projected is _DROP:
                omitted.append(child_path)
                continue
            result.append(projected)
        return result
    if isinstance(value, str) and _looks_like_local_path(value):
        return _DROP
    return canonical_data(value)


def _forbidden_key(value: str) -> bool:
    normalized = str(value).strip().lower()
    if normalized in _SAFE_SEMANTIC_PATH_KEYS:
        return False
    if normalized in _FORBIDDEN_KEYS:
        return True
    if normalized.endswith(
        ("_path", "_paths", "_directory", "_file", "_files")
    ):
        return True
    if any(fragment in normalized for fragment in _SENSITIVE_KEY_FRAGMENTS):
        return True
    if normalized.endswith(("_sha256", "_digest")):
        return False
    return False


def _looks_like_local_path(value: str) -> bool:
    normalized = str(value).strip()
    return bool(
        normalized.startswith(("/", "~/", "./", "../", "file://"))
        or _ABSOLUTE_WINDOWS_PATH.match(normalized)
        or any(
            marker in normalized
            for marker in ("/Users/", "/home/", "/private/", "/tmp/", "/var/")
        )
    )


def _validate_oracle_projection(value: Any) -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            if _forbidden_key(str(key)):
                raise ContractError("host oracle projection contains a forbidden field")
            _validate_oracle_projection(child)
    elif isinstance(value, (tuple, list)):
        for child in value:
            _validate_oracle_projection(child)
    elif isinstance(value, str) and _looks_like_local_path(value):
        raise ContractError("host oracle projection contains a local path")


__all__ = [
    "HostOracleInputBundleV1",
    "HostToolObservationV1",
    "build_host_oracle_input_bundle",
]
