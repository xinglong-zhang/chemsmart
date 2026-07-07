from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class ErrorSignature:
    kind: str
    line: str
    line_no: int | None
    severity: str = "error"


_PATTERNS: list[tuple[str, re.Pattern[str], str]] = [
    (
        "oom_killed",
        re.compile(
            r"(oom-kill|out of memory|killed process \d+|memory cgroup out "
            r"of memory|detected \d+ oom-kill event)",
            re.IGNORECASE,
        ),
        "error",
    ),
    (
        "walltime_exceeded",
        re.compile(
            r"(due to time limit|walltime.*exceed|time limit exceeded|"
            r"job exceeded.*walltime|reached its time limit|cpu time limit "
            r"exceeded)",
            re.IGNORECASE,
        ),
        "error",
    ),
    (
        "missing_module",
        re.compile(
            r"(module:?\s+command not found|unable to locate a modulefile|"
            r"modulefile.*not found|no module named )",
            re.IGNORECASE,
        ),
        "error",
    ),
    (
        "node_failure",
        re.compile(
            r"(node failure|launch failed requeued held|lost communication "
            r"with node|exec_host.*down|node .* not responding|hardware "
            r"failure)",
            re.IGNORECASE,
        ),
        "error",
    ),
    (
        "scheduler_reject",
        re.compile(
            r"(batch job submission failed|job rejected|job exceeds queue "
            r"resource limits|invalid account|unable to run job|denied by "
            r"policy|submission rejected)",
            re.IGNORECASE,
        ),
        "error",
    ),
    (
        "segfault",
        re.compile(
            r"(segmentation fault|sigsegv|segfault at )",
            re.IGNORECASE,
        ),
        "error",
    ),
]


def summarize_log(
    text: str,
    max_signatures: int = 50,
) -> list[ErrorSignature]:
    """Summarize common HPC failure signatures from log text."""

    if max_signatures <= 0:
        return []

    signatures: list[ErrorSignature] = []
    last_seen_by_kind: dict[str, int] = {}

    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        for kind, pattern, severity in _PATTERNS:
            if pattern.search(raw_line) is None:
                continue

            last_line_no = last_seen_by_kind.get(kind)
            if last_line_no is not None and line_no - last_line_no <= 3:
                break

            signatures.append(
                ErrorSignature(
                    kind=kind,
                    line=raw_line.strip(),
                    line_no=line_no,
                    severity=severity,
                )
            )
            last_seen_by_kind[kind] = line_no
            if len(signatures) >= max_signatures:
                return signatures
            break

    return signatures
