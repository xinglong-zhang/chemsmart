from __future__ import annotations

import re
from collections import Counter
from typing import Any

from chemsmart.agent.harness.models import InvariantIssue, InvariantResult

RULE_ID = "gaussian.ts.route"
_OPT_BLOCK_RE = re.compile(r"\bopt\s*=\s*\((?P<body>[^)]*)\)", re.IGNORECASE)


def check_gaussian_ts_route(
    route: str | None,
    *,
    inputfile: str | None = None,
    result_index: int | None = None,
) -> InvariantResult:
    issues: list[InvariantIssue] = []
    evidence: dict[str, Any] = {
        "inputfile": inputfile,
        "result_index": result_index,
        "route": route,
    }
    if route is None:
        issues.append(_reject("Gaussian TS route line is missing", evidence))
        return _result(issues, evidence)

    if "[" in route or "]" in route:
        issues.append(
            _reject(
                "Gaussian TS route contains a Python list/string leak",
                evidence,
            )
        )

    match = _OPT_BLOCK_RE.search(route)
    if match is None:
        issues.append(
            _reject("Gaussian TS route missing opt=(...) block", evidence)
        )
        return _result(issues, evidence)

    tokens = _opt_tokens(match.group("body"))
    token_counts = Counter(tokens)
    evidence = {**evidence, "opt_tokens": tokens}

    if "ts" not in token_counts:
        issues.append(_reject("Gaussian TS route missing ts option", evidence))

    duplicate_tokens = [
        token for token, count in token_counts.items() if count > 1
    ]
    if duplicate_tokens:
        issues.append(
            _reject(
                "Gaussian TS route contains duplicate opt keyword(s)",
                {**evidence, "duplicate_tokens": duplicate_tokens},
            )
        )

    if "calcfc" in token_counts and "calcall" in token_counts:
        issues.append(
            _reject(
                "Gaussian TS route must not combine calcfc and calcall",
                evidence,
            )
        )

    if "noeigentest" not in token_counts:
        issues.append(
            _reject("Gaussian TS route missing noeigentest option", evidence)
        )

    return _result(issues, evidence)


def _opt_tokens(body: str) -> list[str]:
    return [
        token.strip().lower()
        for token in re.split(r"[,\s]+", body)
        if token.strip()
    ]


def _reject(message: str, evidence: dict[str, Any]) -> InvariantIssue:
    return InvariantIssue(
        rule_id=RULE_ID,
        severity="reject",
        message=message,
        evidence=evidence,
    )


def _result(
    issues: list[InvariantIssue], context: dict[str, Any]
) -> InvariantResult:
    return InvariantResult(
        rule_id=RULE_ID,
        verdict="reject" if issues else "ok",
        issues=tuple(issues),
        context=context,
    )
