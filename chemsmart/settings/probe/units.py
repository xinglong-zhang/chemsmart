"""Normalize scheduler-reported quantities to one canonical unit each.

Every memory value becomes integer kilobytes and every wall time becomes
integer seconds the moment it is parsed (the aiida-core rule); MEM_GB and
NUM_HOURS exist only at YAML-write time. Schedulers disagree wildly about
formats -- SLURM's sinfo prints megabytes as bare numbers with an optional
trailing ``+`` on heterogeneous partitions, PBS prints ``60gb``/``123456kb``,
limits appear as ``UNLIMITED``/``infinite``/``n/a`` -- so both parsers accept
the full observed zoo and refuse loudly on anything else.
"""

from __future__ import annotations

import re

_MEMORY = re.compile(
    r"^\s*(?P<number>\d+(?:\.\d+)?)\s*(?P<unit>[kmgt])?b?\+?\s*$",
    re.IGNORECASE,
)
_NO_LIMIT = {"unlimited", "infinite", "infinity", "none", "n/a", "-", ""}
_KB_PER_UNIT = {"k": 1, "m": 1024, "g": 1024**2, "t": 1024**3}


class ProbeUnitError(ValueError):
    """A scheduler-reported quantity could not be understood."""


def parse_memory_kb(text: object, *, default_unit: str = "k") -> int | None:
    """Memory string -> integer kB; ``None`` for an explicit no-limit.

    ``default_unit`` names the unit of a bare number: SLURM's ``sinfo %m``
    prints megabytes (``default_unit="m"``); PBS byte counts default to kB.
    """

    value = str(text).strip().lower()
    if value in _NO_LIMIT:
        return None
    match = _MEMORY.match(value)
    if match is None:
        raise ProbeUnitError(f"unrecognized memory value: {text!r}")
    unit = (match.group("unit") or default_unit).lower()
    if unit not in _KB_PER_UNIT:
        raise ProbeUnitError(f"unrecognized memory unit in: {text!r}")
    return int(float(match.group("number")) * _KB_PER_UNIT[unit])


_TIME_FORMS = (
    # days-hours[:minutes[:seconds]]  (SLURM "2-00:00:00", "1-12")
    re.compile(
        r"^(?P<days>\d+)-(?P<hours>\d+)(?::(?P<minutes>\d+)"
        r"(?::(?P<seconds>\d+))?)?$"
    ),
    # hours:minutes:seconds  (SLURM and PBS "48:00:00")
    re.compile(r"^(?P<hours>\d+):(?P<minutes>\d+):(?P<seconds>\d+)$"),
    # minutes[:seconds]  (SLURM "30", "30:15")
    re.compile(r"^(?P<minutes>\d+)(?::(?P<seconds>\d+))?$"),
)


def parse_walltime_seconds(text: object) -> int | None:
    """Wall-time string -> integer seconds; ``None`` for no stated limit."""

    value = str(text).strip().lower()
    if value in _NO_LIMIT or value == "not_set":
        return None
    for form in _TIME_FORMS:
        match = form.match(value)
        if match is None:
            continue
        parts = {
            name: int(found)
            for name, found in match.groupdict(default="0").items()
        }
        return (
            parts.get("days", 0) * 86400
            + parts.get("hours", 0) * 3600
            + parts.get("minutes", 0) * 60
            + parts.get("seconds", 0)
        )
    raise ProbeUnitError(f"unrecognized wall-time value: {text!r}")


__all__ = ["ProbeUnitError", "parse_memory_kb", "parse_walltime_seconds"]
