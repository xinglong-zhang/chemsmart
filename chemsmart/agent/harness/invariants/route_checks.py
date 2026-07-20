"""Generated-input route invariants for the kinds beyond gaussian.ts.

Each check inspects the *generated* Gaussian/ORCA route line (extracted from the
dry-run input) and returns an :class:`InvariantResult`. Grounded in the real
route builders:

- ORCA TS  -> ``! OptTS ...`` (jobs/orca/settings.py); ``ScanTS`` for scants.
- Frequency-> Gaussian route carries a ` freq` token; ORCA route carries ``Freq``
  (jobs/gaussian/settings.py, jobs/orca/settings.py). Exactly one is expected —
  a missing one means the freq calc silently did not happen (the ORCA-freq bug),
  a duplicate one is the historical "freq twice" bug.
- IRC     -> Gaussian route contains ``irc=(...)``; ORCA route contains ``IRC``.
"""

from __future__ import annotations

import re
from typing import Any

from chemsmart.agent.harness.models import InvariantIssue, InvariantResult

_FREQ_RE = re.compile(r"\b(?:num|an)?freq\b", re.IGNORECASE)
_ORCA_TS_RE = re.compile(r"\bopt\s*-?\s*ts\b", re.IGNORECASE)
_ORCA_SCANTS_RE = re.compile(r"\bscan\s*-?\s*ts\b", re.IGNORECASE)
_GAUSSIAN_IRC_RE = re.compile(r"\birc\s*=", re.IGNORECASE)
_ORCA_IRC_RE = re.compile(r"\bIRC\b", re.IGNORECASE)


def _issue(
    rule_id: str, severity: str, message: str, evidence: dict[str, Any]
) -> InvariantIssue:
    return InvariantIssue(
        rule_id=rule_id, severity=severity, message=message, evidence=evidence
    )


def _result(
    rule_id: str, issues: list[InvariantIssue], context: dict[str, Any]
) -> InvariantResult:
    verdict = (
        "reject"
        if any(i.severity == "reject" for i in issues)
        else ("warn" if issues else "ok")
    )
    return InvariantResult(
        rule_id=rule_id, verdict=verdict, issues=tuple(issues), context=context
    )


def check_orca_ts_route(
    route: str | None,
    *,
    inputfile: str | None = None,
    result_index: int | None = None,
) -> InvariantResult:
    rule_id = "orca.ts.route"
    ev = {"inputfile": inputfile, "result_index": result_index, "route": route}
    issues: list[InvariantIssue] = []
    if route is None:
        issues.append(
            _issue(rule_id, "reject", "ORCA TS route line is missing", ev)
        )
        return _result(rule_id, issues, ev)
    n_optts = len(_ORCA_TS_RE.findall(route))
    n_scants = len(_ORCA_SCANTS_RE.findall(route))
    if n_optts == 0 and n_scants == 0:
        issues.append(
            _issue(
                rule_id,
                "reject",
                "ORCA TS route missing OptTS/ScanTS keyword",
                ev,
            )
        )
    if n_optts > 1 or n_scants > 1:
        issues.append(
            _issue(
                rule_id,
                "reject",
                "ORCA TS route contains a duplicate TS keyword",
                {**ev, "optts_count": n_optts, "scants_count": n_scants},
            )
        )
    return _result(rule_id, issues, ev)


def check_freq_route(
    route: str | None,
    *,
    software: str,
    inputfile: str | None = None,
    result_index: int | None = None,
) -> InvariantResult:
    """The job is a frequency job (dispatched only for ``*.freq`` kinds): the
    generated route MUST contain exactly one freq keyword."""
    rule_id = f"{software}.freq.route"
    ev = {
        "inputfile": inputfile,
        "result_index": result_index,
        "route": route,
        "software": software,
    }
    issues: list[InvariantIssue] = []
    if route is None:
        issues.append(
            _issue(
                rule_id,
                "reject",
                "route line is missing for a frequency job",
                ev,
            )
        )
        return _result(rule_id, issues, ev)
    n = len(_FREQ_RE.findall(route))
    if n == 0:
        issues.append(
            _issue(
                rule_id,
                "reject",
                "frequency job but the generated route has no Freq keyword "
                "(it will silently run without frequencies)",
                ev,
            )
        )
    elif n > 1:
        issues.append(
            _issue(
                rule_id,
                "reject",
                "generated route contains a duplicate freq keyword",
                {**ev, "freq_count": n},
            )
        )
    return _result(rule_id, issues, ev)


def check_irc_route(
    route: str | None,
    *,
    software: str,
    inputfile: str | None = None,
    result_index: int | None = None,
) -> InvariantResult:
    rule_id = f"{software}.irc.route"
    ev = {
        "inputfile": inputfile,
        "result_index": result_index,
        "route": route,
        "software": software,
    }
    issues: list[InvariantIssue] = []
    if route is None:
        issues.append(
            _issue(
                rule_id, "reject", "route line is missing for an IRC job", ev
            )
        )
        return _result(rule_id, issues, ev)
    pattern = _ORCA_IRC_RE if software == "orca" else _GAUSSIAN_IRC_RE
    if pattern.search(route) is None:
        issues.append(
            _issue(
                rule_id,
                "reject",
                f"IRC job but the generated {software} route is missing the required "
                f"IRC keyword",
                ev,
            )
        )
    return _result(rule_id, issues, ev)
