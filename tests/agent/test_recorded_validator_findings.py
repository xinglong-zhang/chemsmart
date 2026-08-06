"""A validation failure must record what it objected to, not only how many.

Observed live: an n-butane session was told twice that a correct ORCA project
produced invalid generated input. The durable record held
``critical_finding_count: 2`` and ``status: invalid`` and nothing else, so the
defect -- ``ri_approximation`` being writable into the route but unreadable
back out of it -- had to be reproduced by hand before it could be found.

This is the same gap that ``tool_failed`` had for rejected arguments: a count
or a digest is not reviewable. These tests pin that findings are recorded, and
that recording them cannot turn an event log into an echo channel.
"""

from dataclasses import dataclass
from typing import Any

from chemsmart.agent.tool_runtime import (
    _FINDING_VALUE_CHARS,
    _RECORDED_FINDINGS,
    _public_validator_findings,
)


@dataclass(frozen=True)
class _Finding:
    rule_id: str
    field: str
    expected: Any
    observed: Any


@dataclass(frozen=True)
class _Validator:
    findings: tuple


def test_the_field_that_disagreed_is_recorded_with_both_values():
    """The exact shape of the live ORCA failure."""

    recorded = _public_validator_findings(
        _Validator(
            (
                _Finding(
                    rule_id="preview.settings.mismatch",
                    field="ri_approximation",
                    expected="none",
                    observed=None,
                ),
            )
        )
    )
    assert len(recorded) == 1
    assert recorded[0]["field"] == "ri_approximation"
    assert recorded[0]["expected"] == "'none'"
    assert recorded[0]["observed"] == "None"
    assert recorded[0]["rule_id"] == "preview.settings.mismatch"


def test_a_long_value_is_truncated_rather_than_echoed():
    recorded = _public_validator_findings(
        _Validator(
            (_Finding("r", "route", "x" * 5000, "y" * 5000),)
        )
    )
    assert len(recorded[0]["expected"]) <= _FINDING_VALUE_CHARS
    assert recorded[0]["expected"].endswith("...")


def test_a_structured_value_is_summarized_by_shape():
    recorded = _public_validator_findings(
        _Validator((_Finding("r", "geometry", [1, 2, 3], {"a": 1}),))
    )
    assert recorded[0]["expected"] == "a list of 3 entries"
    assert recorded[0]["observed"] == "a dict of 1 entries"


def test_many_findings_are_capped_and_the_cap_is_stated():
    recorded = _public_validator_findings(
        _Validator(
            tuple(
                _Finding(f"r{index}", f"f{index}", index, index + 1)
                for index in range(_RECORDED_FINDINGS + 5)
            )
        )
    )
    assert len(recorded) == _RECORDED_FINDINGS + 1
    marker = recorded[-1]
    assert marker["rule_id"] == "record.truncated"
    assert str(_RECORDED_FINDINGS + 5) in marker["expected"]


def test_a_clean_validation_records_nothing():
    assert _public_validator_findings(_Validator(())) == ()


def test_a_validator_without_findings_does_not_raise():
    """Adapters differ; the recorder must not assume a field exists."""

    class _Bare:
        pass

    assert _public_validator_findings(_Bare()) == ()
