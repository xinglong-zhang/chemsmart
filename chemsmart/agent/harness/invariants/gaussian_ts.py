from __future__ import annotations

import re
from collections import Counter
from typing import Any

from chemsmart.agent.harness.models import InvariantIssue, InvariantResult

RULE_ID = "gaussian.ts.route"
_OPT_BLOCK_RE = re.compile(r"\bopt\s*=\s*\((?P<body>[^)]*)\)", re.IGNORECASE)

# Runtime-derived canonical TS triple.
_CANON = {"ts", "calcfc", "noeigentest"}

# Genuine user extras allowed in the TS opt block. Everything else is treated
# as a model leak or malformed route token.
_ALLOWED_EXTRA = re.compile(
    r"^(calcall|cartesian|tight|"
    r"(maxstep|maxcycles|maxcycle|recalcfc|maxmicroiterations)=\d+)$",
    re.IGNORECASE,
)


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

    # A ts job renders exactly one runtime-owned opt=(...) block. A second one
    # (or bare ts/calcfc/noeigentest tokens outside the block) is a model leak
    # via --additional-route-parameters that would corrupt the route.
    opt_blocks = _OPT_BLOCK_RE.findall(route)
    opt_keyword_count = len(re.findall(r"\bopt\b", route, re.IGNORECASE))
    if opt_keyword_count > 1:
        issues.append(
            _reject(
                "Gaussian TS route contains more than one Opt route form",
                {**evidence, "opt_keyword_count": opt_keyword_count},
            )
        )
    if len(opt_blocks) > 1:
        issues.append(
            _reject(
                "Gaussian TS route contains more than one opt=(...) block",
                {**evidence, "opt_block_count": len(opt_blocks)},
            )
        )
    leaked_outside = _leaked_ts_tokens_outside_opt_block(route)
    if leaked_outside:
        issues.append(
            _reject(
                "Gaussian TS route leaks runtime-owned TS token(s) outside "
                "the opt=(...) block",
                {**evidence, "leaked_tokens": sorted(leaked_outside)},
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

    bad_extras = [
        token
        for token in token_counts
        if token not in _CANON and not _ALLOWED_EXTRA.match(token)
    ]
    if bad_extras:
        issues.append(
            _reject(
                "Gaussian TS route contains non-allowlisted extra opt token(s)",
                {**evidence, "bad_extras": sorted(bad_extras)},
            )
        )

    return _result(issues, evidence)


def _opt_tokens(body: str) -> list[str]:
    return [
        token.strip().lower()
        for token in re.split(r"[,\s]+", body)
        if token.strip()
    ]


def _leaked_ts_tokens_outside_opt_block(route: str) -> set[str]:
    """Return runtime-owned TS tokens that appear outside any opt=(...) block."""
    remainder = _OPT_BLOCK_RE.sub(" ", route)
    tokens = {
        token.strip().lower()
        for token in re.split(r"[,\s]+", remainder)
        if token.strip()
    }
    return tokens & _CANON


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
