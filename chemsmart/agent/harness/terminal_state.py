"""Deterministic terminal-state assertions for agent job workflows.

The command semantic gate proves that a command can be parsed and that a safe
input can be rendered.  Submission workflows need one more receipt: the
requested terminal action must have happened and its artifacts must exist.
This module keeps that receipt small, JSON-serializable, and independent of a
provider or hidden reasoning trace.
"""

from __future__ import annotations

from typing import Any, Literal

JsonDict = dict[str, Any]
TerminalStatus = Literal["passed", "failed", "unknown"]

TERMINAL_STATE_SCHEMA_VERSION = 1


def assertion(
    assertion_id: str,
    *,
    expected: Any,
    observed: Any,
    evidence: JsonDict | None = None,
) -> JsonDict:
    """Create one observable assertion row."""

    return {
        "id": str(assertion_id),
        "expected": expected,
        "observed": observed,
        "status": "pass" if expected == observed else "fail",
        "evidence": evidence or {},
    }


def build_terminal_state(
    *,
    action: str,
    command: str,
    assertions: list[JsonDict],
    returncode: int | None = None,
    expected_returncode: int | None = None,
    server: str | None = None,
    scheduler: str | None = None,
    execution_mode: str | None = None,
    safe_execution_mode: str | None = None,
    required_assertion_ids: list[str] | tuple[str, ...] | None = None,
    artifacts: list[JsonDict] | None = None,
) -> JsonDict:
    """Build a terminal receipt from deterministic observations."""

    rows = [dict(row) for row in assertions if isinstance(row, dict)]
    required = tuple(dict.fromkeys(required_assertion_ids or ()))
    observed_ids = {str(row.get("id") or "") for row in rows}
    all_passed = (
        bool(rows)
        and all(row.get("status") == "pass" for row in rows)
        and all(assertion_id in observed_ids for assertion_id in required)
    )
    if expected_returncode is not None and returncode != expected_returncode:
        all_passed = False
    return {
        "schema_version": TERMINAL_STATE_SCHEMA_VERSION,
        "action": action,
        "command": command,
        "status": "passed" if all_passed else "failed",
        "all_passed": all_passed,
        "returncode": returncode,
        "expected_returncode": expected_returncode,
        "server": server,
        "scheduler": scheduler,
        "execution_mode": execution_mode,
        "safe_execution_mode": safe_execution_mode,
        "required_assertion_ids": list(required),
        "assertions": rows,
        "artifacts": [
            dict(item) for item in (artifacts or []) if isinstance(item, dict)
        ],
    }


def validate_terminal_state(state: Any) -> list[str]:
    """Return stable rule IDs for an invalid or incomplete terminal receipt."""

    if not isinstance(state, dict):
        return ["terminal_state.missing"]
    issues = _header_issues(state)
    assertions = state.get("assertions")
    if not isinstance(assertions, list) or not assertions:
        issues.append("terminal_state.assertions_missing")
        assertions = []
    seen: set[str] = set()
    issues.extend(_assertion_row_issues(assertions, seen))
    required = state.get("required_assertion_ids") or []
    issues.extend(_required_assertion_issues(required, seen))
    returncode = state.get("returncode")
    expected_returncode = state.get("expected_returncode")
    issues.extend(_returncode_issues(returncode, expected_returncode))
    computed = (
        bool(assertions)
        and all(
            isinstance(row, dict) and row.get("status") == "pass"
            for row in assertions
        )
        and not (
            expected_returncode is not None
            and returncode != expected_returncode
        )
        and all(assertion_id in seen for assertion_id in required)
    )
    if state.get("all_passed") is not computed:
        issues.append("terminal_state.aggregate_mismatch")
    if state.get("status") != ("passed" if computed else "failed"):
        issues.append("terminal_state.status_mismatch")
    return sorted(set(issues))


def _header_issues(state: JsonDict) -> list[str]:
    issues: list[str] = []
    if state.get("schema_version") != TERMINAL_STATE_SCHEMA_VERSION:
        issues.append("terminal_state.schema_version")
    if not str(state.get("action") or "").strip():
        issues.append("terminal_state.action_missing")
    if not str(state.get("command") or "").strip():
        issues.append("terminal_state.command_missing")
    return issues


def _assertion_row_issues(assertions: list[Any], seen: set[str]) -> list[str]:
    issues: list[str] = []
    for row in assertions:
        if not isinstance(row, dict):
            issues.append("terminal_state.assertion_malformed")
            continue
        assertion_id = str(row.get("id") or "")
        if not assertion_id:
            issues.append("terminal_state.assertion_id_missing")
        elif assertion_id in seen:
            issues.append("terminal_state.assertion_id_duplicate")
        seen.add(assertion_id)
        if row.get("status") not in {"pass", "fail"}:
            issues.append("terminal_state.assertion_status_invalid")
        if row.get("status") == "fail":
            issues.append(f"terminal_state.failed:{assertion_id or 'unknown'}")
    return issues


def _required_assertion_issues(required: Any, seen: set[str]) -> list[str]:
    if not isinstance(required, list) or any(
        not isinstance(item, str) or not item for item in required
    ):
        return ["terminal_state.required_assertions_invalid"]
    return [
        f"terminal_state.required_missing:{assertion_id}"
        for assertion_id in required
        if assertion_id not in seen
    ]


def _returncode_issues(returncode: Any, expected_returncode: Any) -> list[str]:
    issues: list[str] = []
    if expected_returncode is not None and not isinstance(
        expected_returncode, int
    ):
        issues.append("terminal_state.expected_returncode_invalid")
    if (
        isinstance(expected_returncode, int)
        and isinstance(returncode, int)
        and returncode != expected_returncode
    ):
        issues.append("terminal_state.expected_returncode_mismatch")
    if (
        expected_returncode is None
        and isinstance(returncode, int)
        and returncode != 0
    ):
        issues.append("terminal_state.returncode_nonzero")
    return issues


def terminal_state_is_positive(state: Any) -> bool:
    """Whether a receipt is complete and all assertions passed."""

    return (
        isinstance(state, dict)
        and not validate_terminal_state(state)
        and state.get("status") == "passed"
        and state.get("all_passed") is True
        and state.get("returncode") == 0
        and state.get("expected_returncode") in {None, 0}
    )


__all__ = [
    "TERMINAL_STATE_SCHEMA_VERSION",
    "assertion",
    "build_terminal_state",
    "terminal_state_is_positive",
    "validate_terminal_state",
]
