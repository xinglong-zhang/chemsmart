"""Deterministic intent-preservation checks for submitted CLI commands.

The semantic gate answers whether a command is executable.  This module checks
the separate question required by training-data quality: whether the command
preserves the research request that produced it.  A command can be syntactically
valid and successfully submitted while silently changing the input molecule,
program, job type, state, or server.
"""

from __future__ import annotations

import re
from pathlib import PurePath
from typing import Any

from chemsmart.agent.harness.intent import ObservedIntent
from chemsmart.agent.harness.value_equivalence import structured_sequence
from chemsmart.agent.model_command_parser import parse_model_command

JsonDict = dict[str, Any]


def _text(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _same_path(observed: str | None, expected: str | None) -> bool:
    if observed is None or expected is None:
        return observed == expected
    return str(PurePath(observed)) == str(PurePath(expected))


def _option_value(tokens: list[str], names: set[str]) -> str | None:
    """Read an option value that is outside the generic command parser."""

    for index, token in enumerate(tokens):
        for name in names:
            if token == name and index + 1 < len(tokens):
                return tokens[index + 1]
            prefix = f"{name}="
            if token.startswith(prefix):
                return token[len(prefix) :]
    return None


def _equivalent_value(expected: Any, observed: Any) -> bool:
    if expected is None or observed is None:
        return expected == observed
    expected_text = str(expected).strip()
    observed_text = str(observed).strip()
    expected_sequence = structured_sequence(expected)
    if expected_sequence is not None:
        return expected_sequence == structured_sequence(observed)
    try:
        return float(expected_text) == float(observed_text)
    except ValueError:
        expected_compact = re.sub(r"[\s\[\](){}]", "", expected_text)
        expected_numbers = re.findall(r"-?\d+(?:\.\d+)?", expected_text)
        observed_numbers = re.findall(r"-?\d+(?:\.\d+)?", observed_text)
        if expected_numbers and observed_numbers:
            return expected_numbers == observed_numbers
        observed_compact = re.sub(r"[\s\[\](){}]", "", observed_text)
        return expected_compact == observed_compact


def _append_assertion(
    rows: list[JsonDict],
    command: str,
    assertion_id: str,
    expected_value: Any,
    observed_value: Any,
) -> None:
    rows.append(
        {
            "id": assertion_id,
            "expected": expected_value,
            "observed": observed_value,
            "status": (
                "pass"
                if _equivalent_value(expected_value, observed_value)
                else "fail"
            ),
            "evidence": {"command": command},
        }
    )


def _append_core_assertions(
    rows: list[JsonDict],
    command: str,
    expected: JsonDict,
    parsed: Any,
    observed_intent: ObservedIntent,
) -> None:
    _append_assertion(rows, command, "sub.intent_action", "sub", parsed.action)
    if "kind" in expected:
        _append_assertion(
            rows,
            command,
            "sub.intent_kind",
            _text(expected["kind"]),
            observed_intent.kind,
        )
    if "project" in expected:
        _append_assertion(
            rows,
            command,
            "sub.intent_project",
            _text(expected["project"]),
            parsed.project,
        )

    for field, assertion_id in (
        ("program", "sub.intent_program"),
        ("job", "sub.intent_job"),
        ("server", "sub.intent_server"),
        ("charge", "sub.intent_charge"),
        ("multiplicity", "sub.intent_multiplicity"),
    ):
        if field not in expected:
            continue
        observed = _observed_core_value(parsed, field)
        _append_assertion(
            rows,
            command,
            assertion_id,
            _text(expected[field]),
            observed,
        )


def _observed_core_value(parsed: Any, field: str) -> str | None:
    observed = _text(getattr(parsed, field))
    if observed is not None or "qmmm" not in parsed.tokens:
        return observed
    if field == "charge":
        return _option_value(parsed.tokens, {"-ct", "--charge-total"})
    if field == "multiplicity":
        return _option_value(parsed.tokens, {"-mt", "--mult-total"})
    return observed


def _append_filename_assertion(
    rows: list[JsonDict], command: str, expected: JsonDict, parsed: Any
) -> None:
    if "filename" not in expected:
        return
    expected_filename = _text(expected["filename"])
    observed_filename = _text(parsed.filename)
    _append_assertion(
        rows,
        command,
        "sub.intent_filename",
        expected_filename,
        observed_filename
        if _same_path(observed_filename, expected_filename)
        else observed_filename,
    )


def _append_route_assertion(
    rows: list[JsonDict], command: str, expected: JsonDict, parsed: Any
) -> None:
    if "route_contains" not in expected:
        return
    route = (_text(parsed.route_parameters) or "").lower()
    required = (_text(expected["route_contains"]) or "").lower()
    rows.append(
        {
            "id": "sub.intent_route_contains",
            "expected": required,
            "observed": route,
            "status": "pass" if required in route else "fail",
            "evidence": {"command": command},
        }
    )


def _append_structural_assertions(
    rows: list[JsonDict], command: str, expected: JsonDict, parsed: Any
) -> None:
    for field, assertion_id in (
        ("ending_xyzfile", "sub.intent_ending_xyzfile"),
        ("nimages", "sub.intent_nimages"),
        ("coordinates", "sub.intent_coordinates"),
        ("num_steps", "sub.intent_num_steps"),
        ("step_size", "sub.intent_step_size"),
        ("dist_start", "sub.intent_dist_start"),
        ("dist_end", "sub.intent_dist_end"),
        ("constrained_coordinates", "sub.intent_constrained_coordinates"),
    ):
        if field not in expected:
            continue
        observed = _structural_value(parsed, field)
        expected_value = expected[field]
        observed_value = observed
        if field not in {"coordinates", "constrained_coordinates"}:
            expected_value = _text(expected_value)
            observed_value = _text(observed_value)
        _append_assertion(
            rows,
            command,
            assertion_id,
            expected_value,
            observed_value,
        )


def _structural_value(parsed: Any, field: str) -> Any:
    observed = parsed.structural_options.get(field)
    if observed is None and field == "num_steps":
        return parsed.structural_options.get("num_steps_or_every_n_points")
    if observed is None and field == "step_size":
        return parsed.structural_options.get("step_size_or_solv")
    return observed


def _append_required_token_assertions(
    rows: list[JsonDict], command: str, expected: JsonDict, parsed: Any
) -> None:
    if "required_tokens" not in expected:
        return
    tokens = set(parsed.tokens)
    for token in expected["required_tokens"] or []:
        token_text = str(token)
        _append_assertion(
            rows,
            command,
            f"sub.intent_required_token:{token_text}",
            True,
            token_text in tokens,
        )


def build_sub_intent_assertions(
    command: str,
    expected: JsonDict,
    *,
    cwd: str | None = None,
) -> list[JsonDict]:
    """Return observable assertions for a ``chemsmart sub`` intent.

    ``expected`` may contain ``program``, ``job``, ``server``, ``filename``,
    ``charge``, ``multiplicity``, ``route_contains``, ``ending_xyzfile``,
    ``nimages``, and ``required_tokens``.  Missing expected keys are not
    asserted, so callers can use the same helper for simple and specialized
    chemistry workflows.
    """

    parsed = parse_model_command(command, cwd=cwd)
    observed_intent = ObservedIntent.from_command(command, cwd=cwd)
    rows: list[JsonDict] = []
    _append_core_assertions(rows, command, expected, parsed, observed_intent)
    _append_filename_assertion(rows, command, expected, parsed)
    _append_route_assertion(rows, command, expected, parsed)
    _append_structural_assertions(rows, command, expected, parsed)
    _append_required_token_assertions(rows, command, expected, parsed)
    return rows


__all__ = ["build_sub_intent_assertions"]
