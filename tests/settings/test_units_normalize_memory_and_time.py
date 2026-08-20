"""One canonical unit each: kilobytes and seconds, parsed at the boundary.

Schedulers disagree about formats; the parsers accept the observed zoo and
refuse loudly on anything else, so a mis-read limit can never silently
become a wrong resource request.
"""

from __future__ import annotations

import pytest

from chemsmart.settings.probe import (
    ProbeUnitError,
    parse_memory_kb,
    parse_walltime_seconds,
)


@pytest.mark.parametrize(
    ("text", "default_unit", "expected_kb"),
    (
        ("60000", "m", 60000 * 1024),  # sinfo %m prints megabytes bare
        ("60000+", "m", 60000 * 1024),  # heterogeneous partition marker
        ("375G", "k", 375 * 1024 * 1024),
        ("1.5T", "k", int(1.5 * 1024**3)),
        ("187000M", "k", 187000 * 1024),
        ("123456kb", "k", 123456),  # PBS spelling
        ("60gb", "k", 60 * 1024 * 1024),
        ("2048", "k", 2048),
    ),
)
def test_memory_forms_normalize_to_kb(text, default_unit, expected_kb):
    assert parse_memory_kb(text, default_unit=default_unit) == expected_kb


@pytest.mark.parametrize("text", ("UNLIMITED", "infinite", "n/a", "-", ""))
def test_no_limit_memory_is_none_not_zero(text):
    assert parse_memory_kb(text) is None


@pytest.mark.parametrize("text", ("sixty", "12q", "1..5G", "G"))
def test_nonsense_memory_refuses_loudly(text):
    with pytest.raises(ProbeUnitError):
        parse_memory_kb(text)


@pytest.mark.parametrize(
    ("text", "expected_seconds"),
    (
        ("2-00:00:00", 172800),  # SLURM MaxTime on this host
        ("1-12", 129600),
        ("1-2:30", 95400),
        ("48:00:00", 172800),  # PBS walltime
        ("00:30:00", 1800),
        ("30", 1800),  # bare SLURM minutes
        ("30:15", 1815),
    ),
)
def test_walltime_forms_normalize_to_seconds(text, expected_seconds):
    assert parse_walltime_seconds(text) == expected_seconds


@pytest.mark.parametrize("text", ("UNLIMITED", "infinite", "n/a", "NOT_SET"))
def test_no_limit_walltime_is_none(text):
    assert parse_walltime_seconds(text) is None


@pytest.mark.parametrize("text", ("1:2:3:4", "one hour", "1-", "-5"))
def test_nonsense_walltime_refuses_loudly(text):
    with pytest.raises(ProbeUnitError):
        parse_walltime_seconds(text)
