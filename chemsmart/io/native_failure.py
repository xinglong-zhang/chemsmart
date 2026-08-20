"""Bounded, parser-owned summaries of native program failures.

The summaries in this module are deliberately much smaller than a native
output file, and they carry two different kinds of statement.

``diagnostic_lines`` are parser-owned canonical templates keyed by a stable
``error_class``.  They are host statements: short, predictable, and safe to
match on.

``engine_lines`` are the program's *own* words about the failure, quoted
verbatim and bounded.  They exist because the canonical vocabulary is closed:
a failure mode nobody anticipated classifies as ``native_runtime`` and its
template says only that an error occurred, which is no more useful than
silence.  A run once died because canonical CCSD under the RIJK approximation
is not size-consistent -- a qualitative error for exactly the interaction
energy being computed -- and the program said so plainly while the host
reported nothing.  Quoting the engine is not re-implementing its rules, so no
second validation system appears here; the quote is evidence about the run,
never a host claim, and never readiness or validation.

Paths, URLs, and credential-like values are redacted from quoted lines. Input
echoes may appear in a quoted engine line, because what the program objected to
is usually the thing a reader needs to see.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from itertools import chain
from typing import Iterable

_SCHEMA_VERSION = "chemsmart.native-failure-summary.v1"
_MAX_DIAGNOSTICS = 3
_MAX_DIAGNOSTIC_CHARS = 160

#: Programs whose output this module knows how to classify.
_SUPPORTED_PROGRAMS = frozenset({"gaussian", "orca", "xtb"})

#: Engines state the diagnosis and then abort, so the terminating line is
#: frequently the least informative one in the block: ORCA's abort marker
#: literally elides its own reason, and the four lines naming the rule and its
#: remedy come before it.  Quote a bounded window on both sides of the match.
_LEADING_ENGINE_LINES = 4
_TRAILING_ENGINE_LINES = 2

#: Engine lines are quoted rather than templated, so the block needs room for a
#: rule, its consequence and its remedy on either side of the marker.  Still
#: tightly bounded: at most this many lines, each truncated to the shared
#: per-line limit, so the worst case is a little over a kilobyte against a
#: native output of megabytes.
_MAX_ENGINE_LINES = 7


@dataclass(frozen=True)
class NativeFailureSummaryV1:
    """One bounded native-failure observation produced by an output parser."""

    schema_version: str
    program: str
    termination_state: str
    error_class: str
    diagnostic_lines: tuple[str, ...]
    engine_lines: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != _SCHEMA_VERSION:
            raise ValueError("unsupported native failure summary")
        if self.program not in _SUPPORTED_PROGRAMS:
            raise ValueError("unsupported native failure program")
        if self.termination_state not in {"error_termination", "incomplete"}:
            raise ValueError("unsupported native termination state")
        if not re.fullmatch(r"[a-z][a-z0-9_]*", self.error_class):
            raise ValueError("native failure class is not stable")
        for label, lines, limit in (
            ("diagnostics", self.diagnostic_lines, _MAX_DIAGNOSTICS),
            ("engine lines", self.engine_lines, _MAX_ENGINE_LINES),
        ):
            if len(lines) > limit:
                raise ValueError(f"native failure {label} are not bounded")
            if len(set(lines)) != len(lines):
                raise ValueError(f"native failure {label} are not unique")
            for line in lines:
                if not line or "\n" in line or "\r" in line:
                    raise ValueError(f"native failure {label} is not one line")
                if len(line) > _MAX_DIAGNOSTIC_CHARS:
                    raise ValueError(f"native failure {label} is too long")
        for line in self.engine_lines:
            if _REDACTION_REQUIRED.search(line):
                raise ValueError("native failure engine line is not redacted")

    def as_dict(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "program": self.program,
            "termination_state": self.termination_state,
            "error_class": self.error_class,
            "diagnostic_lines": self.diagnostic_lines,
            "engine_lines": self.engine_lines,
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

_XTB_RULES = (
    (
        "scf_convergence",
        (
            re.compile(r"\bSCC\b.*\bnot converged\b", re.I),
            re.compile(r"\bself-consistent charge\b.*\bfailed\b", re.I),
            re.compile(r"\bconvergence\b.*\bnot\b.*\breached\b", re.I),
        ),
    ),
    (
        "geometry_optimization",
        (
            re.compile(r"\bFAILED TO CONVERGE GEOMETRY\b", re.I),
            re.compile(r"\bgeometry optimization\b.*\bfailed\b", re.I),
        ),
    ),
    (
        "input_syntax",
        (
            re.compile(
                r"\bcould not read\b.*\b(?:file|input|geometry)\b", re.I
            ),
            re.compile(r"\bunknown\b.*\b(?:argument|option|keyword)\b", re.I),
        ),
    ),
    (
        "method_unavailable",
        (
            re.compile(
                r"\bparameter(?:s|isation|ization)?\b.*\bnot\b.*\bavailable\b",
                re.I,
            ),
            re.compile(r"\bno parameters\b.*\belement\b", re.I),
        ),
    ),
)

_ORCA_NORMAL = re.compile(r"ORCA TERMINATED NORMALLY", re.I)
_ORCA_ERROR = re.compile(r"ORCA finished by error termination", re.I)
_ORCA_ABORT = re.compile(r"\bError \(ORCA_MAIN\):.*\baborting the run\b", re.I)
_GAUSSIAN_NORMAL = re.compile(r"Normal termination of Gaussian", re.I)
_GAUSSIAN_ERROR = re.compile(r"Error termination via", re.I)
_XTB_NORMAL = re.compile(r"\*\s*finished run on", re.I)
_XTB_ERROR = re.compile(r"\[ERROR\]|\babnormal termination of xtb\b", re.I)

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
    "xtb": {
        "scf_convergence": ("xTB self-consistent charges did not converge.",),
        "geometry_optimization": (
            "xTB geometry optimization did not converge.",
        ),
        "input_syntax": ("xTB rejected its input or command line.",),
        "method_unavailable": (
            "xTB has no parameters for the requested method or element.",
        ),
        "native_runtime": ("xTB reported an error termination.",),
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


def summarize_xtb_native_failure(
    output_lines: Iterable[str],
    *,
    diagnostic_lines: Iterable[str] = (),
) -> NativeFailureSummaryV1 | None:
    """Classify an abnormal xTB run.

    xTB is one of the three programs this release actually executes, so a run
    of it that fails deserves the same account as ORCA or Gaussian.
    """

    return _summarize(
        program="xtb",
        lines=chain(output_lines, diagnostic_lines),
        rules=_XTB_RULES,
        normal_pattern=_XTB_NORMAL,
        error_patterns=(_XTB_ERROR,),
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
    engine_lines: list[str] = []
    recent: list[str] = []
    trailing = 0

    for line_index, raw_line in enumerate(lines):
        line = " ".join(str(raw_line).strip().split())
        if not line:
            continue
        if normal_pattern.search(line):
            last_normal_index = line_index
        matched = False
        if any(pattern.search(line) for pattern in error_patterns):
            explicit_error = True
            last_error_index = line_index
            matched = True
        for error_class, patterns in rules:
            if any(pattern.search(line) for pattern in patterns):
                explicit_error = True
                last_error_index = line_index
                matches[error_class] = True
                matched = True
        # Quote a window around the line the program failed on.  The lines
        # before it carry the diagnosis and usually the remedy; the matched
        # line itself is often only the abort marker.
        if matched:
            for earlier in recent:
                _append_engine_line(engine_lines, earlier)
            _append_engine_line(engine_lines, _redact(line))
            trailing = _TRAILING_ENGINE_LINES
        elif trailing:
            _append_engine_line(engine_lines, _redact(line))
            trailing -= 1
        if _carries_a_diagnosis(line):
            recent.append(_redact(line))
            del recent[:-_LEADING_ENGINE_LINES]

    if last_normal_index >= 0 and last_normal_index > last_error_index:
        return None

    error_class = next(
        (candidate for candidate, _patterns in rules if matches[candidate]),
        "native_runtime" if explicit_error else "incomplete_output",
    )
    diagnostic_candidates = list(_CANONICAL_DIAGNOSTICS[program][error_class])
    diagnostics: list[str] = []
    for line in diagnostic_candidates:
        _append_bounded(diagnostics, line)

    if not engine_lines and error_class == "incomplete_output":
        # An output that simply stops has no marker to quote around, so the
        # window above never opens and the failure reaches the session as a
        # class with no reason attached. The last substantive lines are where
        # the engine got to, and that is usually where it said why.
        #
        # Observed on a real approved run: a relaxed scan ended on "GSTEP:
        # could not impose initial constraints. Please check your input and
        # try again!" -- the engine naming the exact problem -- and the typed
        # evidence carried `engine_lines: []`. A session told only
        # "incomplete_output" cannot repair that; told what ORCA said, it can.
        for line in recent:
            _append_engine_line(engine_lines, line)

    return NativeFailureSummaryV1(
        schema_version=_SCHEMA_VERSION,
        program=program,
        termination_state=(
            "error_termination" if explicit_error else "incomplete"
        ),
        error_class=error_class,
        diagnostic_lines=tuple(diagnostics),
        engine_lines=tuple(engine_lines),
    )


def _append_bounded(values: list[str], line: str) -> None:
    if line and line not in values and len(values) < _MAX_DIAGNOSTICS:
        values.append(line[:_MAX_DIAGNOSTIC_CHARS])


def _append_engine_line(values: list[str], line: str) -> None:
    if line and line not in values and len(values) < _MAX_ENGINE_LINES:
        values.append(line[:_MAX_DIAGNOSTIC_CHARS])


def _carries_a_diagnosis(line: str) -> bool:
    """Whether a preceding line could hold the reason, or is only framing.

    Program output around a failure is full of banner frames and separator
    rules.  They are worth skipping when looking *backwards* for the
    diagnosis, because a fixed window spent on framing crowds out a remedy the
    engine printed after the error instead of before it.
    """

    substance = [character for character in line if character.isalnum()]
    return len(substance) >= 8 and len(substance) >= 0.4 * len(line)


#: A URL, an absolute filesystem path of at least two segments, a home-relative
#: path, an e-mail address, an authorization scheme, or an assignment whose name
#: reads like a credential.  Anything matching this must not survive into a
#: quoted line.
#:
#: A credential assignment is consumed to the end of the line rather than to the
#: next space, because a secret may contain spaces and half of one is still one.
_REDACTION_REQUIRED = re.compile(
    r"""
    https?://\S+
    | (?:/[\w.+-]+){2,}/?
    | ~/[\w.+-/]+
    | \b[\w.+-]+@[\w-]+\.[\w.-]+\b
    | \b(?:bearer|basic)\s+\S+.*
    | \b(?:authorization|api[_-]?key|access[_-]?token|token|secret
        |password|passwd|credential)s?\b\s*[:=]\s*.*
    """,
    re.IGNORECASE | re.VERBOSE,
)


def _redact(line: str) -> str:
    """Remove locators and credential-like values from an engine line.

    What the program objected to is usually the thing a reader needs, so this
    removes only what identifies *this host or account* rather than the
    science: URLs, absolute and home-relative paths, e-mail addresses, and
    assignments whose name reads like a credential.
    """

    return _REDACTION_REQUIRED.sub("<redacted>", line)


__all__ = [
    "NativeFailureSummaryV1",
    "summarize_gaussian_native_failure",
    "summarize_orca_native_failure",
    "summarize_xtb_native_failure",
]
