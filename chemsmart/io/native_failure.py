"""Bounded, parser-owned summaries of native program failures.

The summaries in this module are deliberately much smaller than a native
output file. They retain only a stable failure class and parser-owned canonical
diagnostic templates. Matched raw lines, input echoes, native Link identifiers,
paths, URLs, and credential-like values never enter the public record.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import chain
import re
from typing import Iterable


_SCHEMA_VERSION = "chemsmart.native-failure-summary.v1"
_MAX_DIAGNOSTICS = 3
_MAX_DIAGNOSTIC_CHARS = 160


@dataclass(frozen=True)
class NativeFailureSummaryV1:
    """One bounded native-failure observation produced by an output parser."""

    schema_version: str
    program: str
    termination_state: str
    error_class: str
    diagnostic_lines: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.schema_version != _SCHEMA_VERSION:
            raise ValueError("unsupported native failure summary")
        if self.program not in {"orca", "gaussian"}:
            raise ValueError("unsupported native failure program")
        if self.termination_state not in {"error_termination", "incomplete"}:
            raise ValueError("unsupported native termination state")
        if not re.fullmatch(r"[a-z][a-z0-9_]*", self.error_class):
            raise ValueError("native failure class is not stable")
        if len(self.diagnostic_lines) > _MAX_DIAGNOSTICS:
            raise ValueError("native failure diagnostics are not bounded")
        if len(set(self.diagnostic_lines)) != len(self.diagnostic_lines):
            raise ValueError("native failure diagnostics are not unique")
        for line in self.diagnostic_lines:
            if not line or "\n" in line or "\r" in line:
                raise ValueError("native failure diagnostic is not one line")
            if len(line) > _MAX_DIAGNOSTIC_CHARS:
                raise ValueError("native failure diagnostic is too long")

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "program": self.program,
            "termination_state": self.termination_state,
            "error_class": self.error_class,
            "diagnostic_lines": self.diagnostic_lines,
        }


_ORCA_RULES = (
    (
        "auxiliary_basis",
        (
            re.compile(
                r"not an appropriate auxiliary basis set for correlated methods",
                re.IGNORECASE,
            ),
            re.compile(r"\bauxc\b.*\bbasis\b.*\b(?:missing|required)\b", re.I),
            re.compile(r"\bno\b.*\bauxc\b.*\bbasis\b", re.I),
        ),
    ),
    (
        "mpi_runtime",
        (
            re.compile(r"\bMPI_Type_match_size\b", re.I),
            re.compile(r"\bMPI_ERR_[A-Z_]+\b", re.I),
            re.compile(r"\bPMIX ERROR\b", re.I),
            re.compile(r"\bORTE_ERROR_LOG\b", re.I),
        ),
    ),
    (
        "scf_convergence",
        (
            re.compile(
                r"\bSCF\b.*\b(?:NOT CONVERGED|CONVERGENCE FAILURE)\b", re.I
            ),
            re.compile(r"\bSCF\b.*\bfailed to converge\b", re.I),
        ),
    ),
    (
        "input_syntax",
        (
            re.compile(r"\b(?:unrecognized|unknown)\b.*\bkeyword\b", re.I),
            re.compile(r"\berror\b.*\binput\b.*\bline\b", re.I),
        ),
    ),
    (
        "memory_resource",
        (
            re.compile(r"\b(?:out of memory|not enough memory)\b", re.I),
            re.compile(r"\bcannot allocate memory\b", re.I),
        ),
    ),
    (
        "storage_io",
        (
            re.compile(r"\bno space left on device\b", re.I),
            re.compile(r"\b(?:read|write|i/o) error\b", re.I),
        ),
    ),
)

_GAUSSIAN_RULES = (
    (
        "input_syntax",
        (
            re.compile(r"\bQPErr\b", re.I),
            re.compile(r"syntax error was detected in the input line", re.I),
        ),
    ),
    (
        "scf_convergence",
        (
            re.compile(r"Convergence failure -- run terminated", re.I),
            re.compile(r"\bSCF\b.*\bfailed to converge\b", re.I),
        ),
    ),
    (
        "memory_resource",
        (
            re.compile(r"\bgalloc\b.*\b(?:memory|allocate)\b", re.I),
            re.compile(r"\b(?:out of memory|cannot allocate memory)\b", re.I),
        ),
    ),
    (
        "storage_io",
        (
            re.compile(r"\bErroneous write\b", re.I),
            re.compile(r"\bno space left on device\b", re.I),
            re.compile(r"\b(?:read|write|i/o) error\b", re.I),
        ),
    ),
)

_ORCA_NORMAL = re.compile(r"ORCA TERMINATED NORMALLY", re.I)
_ORCA_ERROR = re.compile(r"ORCA finished by error termination", re.I)
_ORCA_ABORT = re.compile(r"\bError \(ORCA_MAIN\):.*\baborting the run\b", re.I)
_GAUSSIAN_NORMAL = re.compile(r"Normal termination of Gaussian", re.I)
_GAUSSIAN_ERROR = re.compile(r"Error termination via", re.I)

_CANONICAL_DIAGNOSTICS = {
    "orca": {
        "auxiliary_basis": (
            "ORCA rejected the auxiliary basis for a correlated method.",
        ),
        "mpi_runtime": ("ORCA subprocess reported an MPI runtime failure.",),
        "scf_convergence": ("ORCA SCF did not converge.",),
        "input_syntax": ("ORCA rejected native input syntax.",),
        "memory_resource": ("ORCA reported insufficient memory.",),
        "storage_io": ("ORCA reported a storage I/O failure.",),
        "native_runtime": ("ORCA reported an error termination.",),
        "incomplete_output": (),
    },
    "gaussian": {
        "input_syntax": ("Gaussian reported QPErr input syntax failure.",),
        "scf_convergence": ("Gaussian SCF did not converge.",),
        "memory_resource": ("Gaussian reported insufficient memory.",),
        "storage_io": ("Gaussian reported a storage I/O failure.",),
        "native_runtime": ("Gaussian reported an error termination.",),
        "incomplete_output": (),
    },
}


def summarize_orca_native_failure(
    output_lines: Iterable[str],
    *,
    diagnostic_lines: Iterable[str] = (),
) -> NativeFailureSummaryV1 | None:
    """Classify an abnormal ORCA output without returning the native log."""

    return _summarize(
        program="orca",
        lines=chain(output_lines, diagnostic_lines),
        rules=_ORCA_RULES,
        normal_pattern=_ORCA_NORMAL,
        error_patterns=(_ORCA_ERROR, _ORCA_ABORT),
    )


def summarize_gaussian_native_failure(
    output_lines: Iterable[str],
) -> NativeFailureSummaryV1 | None:
    """Classify an abnormal Gaussian output without returning input echoes."""

    return _summarize(
        program="gaussian",
        lines=output_lines,
        rules=_GAUSSIAN_RULES,
        normal_pattern=_GAUSSIAN_NORMAL,
        error_patterns=(_GAUSSIAN_ERROR,),
    )


def _summarize(
    *,
    program: str,
    lines: Iterable[str],
    rules: tuple[tuple[str, tuple[re.Pattern[str], ...]], ...],
    normal_pattern: re.Pattern[str],
    error_patterns: tuple[re.Pattern[str], ...],
) -> NativeFailureSummaryV1 | None:
    matches = {error_class: False for error_class, _patterns in rules}
    last_normal_index = -1
    last_error_index = -1
    explicit_error = False

    for line_index, raw_line in enumerate(lines):
        line = " ".join(str(raw_line).strip().split())
        if not line:
            continue
        if normal_pattern.search(line):
            last_normal_index = line_index
        if any(pattern.search(line) for pattern in error_patterns):
            explicit_error = True
            last_error_index = line_index
        for error_class, patterns in rules:
            if any(pattern.search(line) for pattern in patterns):
                explicit_error = True
                last_error_index = line_index
                matches[error_class] = True

    if last_normal_index >= 0 and last_normal_index > last_error_index:
        return None

    error_class = next(
        (candidate for candidate, _patterns in rules if matches[candidate]),
        "native_runtime" if explicit_error else "incomplete_output",
    )
    diagnostic_candidates = list(
        _CANONICAL_DIAGNOSTICS[program][error_class]
    )
    diagnostics: list[str] = []
    for line in diagnostic_candidates:
        _append_bounded(diagnostics, line)

    return NativeFailureSummaryV1(
        schema_version=_SCHEMA_VERSION,
        program=program,
        termination_state=(
            "error_termination" if explicit_error else "incomplete"
        ),
        error_class=error_class,
        diagnostic_lines=tuple(diagnostics),
    )


def _append_bounded(values: list[str], line: str) -> None:
    if line and line not in values and len(values) < _MAX_DIAGNOSTICS:
        values.append(line[:_MAX_DIAGNOSTIC_CHARS])


__all__ = [
    "NativeFailureSummaryV1",
    "summarize_gaussian_native_failure",
    "summarize_orca_native_failure",
]
