"""BSE-backed basis-set name validation for chemsmart harnesses.

This module intentionally validates *names and intent*, not atomic basis
coefficients. Coefficient/source-of-truth data remains in Basis Set Exchange;
chemsmart keeps a small generated catalog so runtime harness checks are fast
and reproducible offline.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Any, Literal

BasisProgram = Literal["gaussian", "orca"]
BasisIntentVerdict = Literal["ok", "warn", "reject", "ask_user"]
BasisRole = Literal["any", "orbital", "jfit", "jkfit", "rifit", "admmfit"]

_DATA_FILE = "bse_basis_catalog.json"
_TOKEN_SPLIT_RE = re.compile(r"[^a-z0-9*+]+")
_AMBIGUOUS_QUALITY_RE = re.compile(
    r"\b(good|reasonable|standard|large|small|better|best|accurate|cheap|"
    r"fast|diffuse|polarized|polarised|double[-\s]*zeta|triple[-\s]*zeta|"
    r"karlsruhe|dunning|pople)\b",
    re.IGNORECASE,
)
_MAX_SEARCH_LIMIT = 20

_FAMILY_TERMS = {
    "ahlrichs": {"ahlrichs", "karlsruhe", "def2", "def"},
    "dunning": {"dunning", "correlation consistent", "cc", "aug cc"},
    "pople": {
        "pople",
        "split valence",
        "split-valence",
        "six thirty one",
        "6-31",
        "631",
    },
}
_QUALITY_TERMS = {
    "double_zeta": {
        "double zeta",
        "double-zeta",
        "dz",
        "svp",
        "valence double",
    },
    "triple_zeta": {
        "triple zeta",
        "triple-zeta",
        "tz",
        "tzv",
        "tzvp",
        "tee zeta",
        "three zeta",
        "valence triple",
    },
    "quadruple_zeta": {"quadruple zeta", "quadruple-zeta", "qz", "qzvp"},
    "diffuse": {"diffuse", "augmented", "aug", "+", "svpd", "tzvpd"},
    "polarized": {
        "polarized",
        "polarised",
        "polarization",
        "polarisation",
        "star",
        "*",
        "vp",
    },
}
_ROLE_TERMS = {
    "jfit": {"jfit", "j fit", "coulomb fitting", "coulomb-fit"},
    "jkfit": {"jkfit", "jk fit", "exchange fitting", "exchange-fit"},
    "rifit": {"rifit", "ri fit", "ri-fit", "resolution of identity", "aux"},
    "admmfit": {"admm", "admmfit", "admm fit"},
}


@dataclass(frozen=True)
class BasisIntentResult:
    verdict: BasisIntentVerdict
    input_text: str
    program: BasisProgram
    canonical_name: str | None = None
    catalog_key: str | None = None
    message: str = ""
    candidates: tuple[str, ...] = ()
    evidence: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "verdict": self.verdict,
            "input_text": self.input_text,
            "program": self.program,
            "canonical_name": self.canonical_name,
            "catalog_key": self.catalog_key,
            "message": self.message,
            "candidates": list(self.candidates),
            "evidence": self.evidence or {},
        }


@lru_cache(maxsize=1)
def load_basis_catalog() -> dict[str, Any]:
    path = files(__package__).joinpath(_DATA_FILE)
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve_basis_name(
    name: str,
    *,
    program: BasisProgram,
    catalog: dict[str, Any] | None = None,
) -> BasisIntentResult:
    """Resolve an explicit basis-set name against a program-specific catalog."""

    catalog = catalog or load_basis_catalog()
    normalized = normalize_basis_name(name)
    aliases = catalog.get("aliases", {})
    entry_key = aliases.get(normalized)
    if entry_key is None:
        return BasisIntentResult(
            verdict="reject",
            input_text=name,
            program=program,
            message="basis set name is not in the BSE-backed chemsmart catalog",
            candidates=_near_matches(normalized, catalog),
            evidence={"normalized": normalized},
        )

    entry = catalog["basis_sets"][entry_key]
    programs = set(entry.get("programs", []))
    if program not in programs:
        return BasisIntentResult(
            verdict="reject",
            input_text=name,
            program=program,
            canonical_name=entry.get("display_name"),
            catalog_key=entry_key,
            message=f"basis set exists in BSE but is not renderable for {program}",
            evidence={
                "normalized": normalized,
                "programs": sorted(programs),
                "family": entry.get("family"),
                "role": entry.get("role"),
            },
        )

    return BasisIntentResult(
        verdict="ok",
        input_text=name,
        program=program,
        canonical_name=entry.get("display_name"),
        catalog_key=entry_key,
        message="basis set name resolved to a BSE canonical entry",
        evidence={
            "normalized": normalized,
            "family": entry.get("family"),
            "role": entry.get("role"),
            "function_types": entry.get("function_types", []),
            "elements_count": len(entry.get("elements", [])),
        },
    )


def check_basis_intent(
    text: str,
    *,
    program: BasisProgram,
    catalog: dict[str, Any] | None = None,
) -> BasisIntentResult:
    """Classify whether user/model text defines a concrete basis-set name.

    This is deliberately conservative. If a phrase sounds like a basis-set
    family or quality request but does not contain a concrete BSE name, the
    result is ``ask_user`` rather than guessing a basis.
    """

    catalog = catalog or load_basis_catalog()
    names = _basis_mentions(text, catalog)
    if len(names) == 1:
        return resolve_basis_name(names[0], program=program, catalog=catalog)
    if len(names) > 1:
        return BasisIntentResult(
            verdict="warn",
            input_text=text,
            program=program,
            message=(
                "multiple concrete basis names were found; caller must decide "
                "whether this is a mixed-basis request"
            ),
            candidates=tuple(names),
            evidence={"mentions": names},
        )
    if _AMBIGUOUS_QUALITY_RE.search(text or ""):
        return BasisIntentResult(
            verdict="ask_user",
            input_text=text,
            program=program,
            message=(
                "basis intent is qualitative or family-level, not a concrete "
                "basis-set name"
            ),
            candidates=tuple(_family_candidates(text, catalog)),
            evidence={"reason": "qualitative_or_family_basis_intent"},
        )
    return BasisIntentResult(
        verdict="reject",
        input_text=text,
        program=program,
        message="no basis-set intent or concrete basis-set name was detected",
    )


def search_basis_sets(
    query: str,
    *,
    program: BasisProgram = "gaussian",
    limit: int = 8,
    role: BasisRole = "any",
) -> dict[str, Any]:
    """Search BSE-backed basis-set names for a short user phrase.

    This tool is intentionally top-k only. It exists so model providers can
    resolve user phrases such as "Karlsruhe triple zeta diffuse" or "RI fit
    for def2-TZVP" without receiving the full BSE catalog in the prompt.
    """

    catalog = load_basis_catalog()
    limit = max(1, min(int(limit or 8), _MAX_SEARCH_LIMIT))
    query = (query or "").strip()
    normalized_query = normalize_basis_name(_expand_common_spoken_forms(query))
    intent = _query_intent(query)
    requested_role = _infer_role(query, role)
    scored: list[tuple[int, str, list[str]]] = []

    for key, entry in catalog.get("basis_sets", {}).items():
        if program not in set(entry.get("programs", [])):
            continue
        if requested_role != "any" and entry.get("role") != requested_role:
            continue

        score, reasons = _score_entry(
            key=key,
            entry=entry,
            normalized_query=normalized_query,
            query_intent=intent,
            requested_role=requested_role,
        )
        if score > 0:
            scored.append((score, key, reasons))

    ranked = sorted(
        scored,
        key=lambda item: (
            -item[0],
            catalog["basis_sets"][item[1]].get("role") != "orbital",
            len(catalog["basis_sets"][item[1]].get("display_name", "")),
            catalog["basis_sets"][item[1]].get("display_name", ""),
        ),
    )
    candidates = [
        _candidate_payload(catalog["basis_sets"][key], score, reasons)
        for score, key, reasons in ranked[:limit]
    ]
    verdict: BasisIntentVerdict = "ok" if candidates else "reject"
    if candidates and _qualitative_intent_only(intent, normalized_query):
        verdict = "ask_user"
    elif len(candidates) > 1 and ranked[0][0] - ranked[1][0] < 15:
        verdict = "warn"

    return {
        "ok": bool(candidates),
        "verdict": verdict,
        "query": query,
        "program": program,
        "normalized_query": normalized_query,
        "requested_role": requested_role,
        "result_count": len(candidates),
        "limit": limit,
        "truncated": len(ranked) > limit,
        "token_policy": "top_k_only; full catalog is never returned",
        "candidates": candidates,
        "message": _search_message(verdict, candidates),
    }


def normalize_basis_name(name: str) -> str:
    text = (name or "").strip().lower()
    text = text.replace("ζ", "zeta")
    # BSE names are punctuation-sensitive in display form, but harness lookup
    # should tolerate user/model spelling variants such as def2 TZVP.
    return "".join(part for part in _TOKEN_SPLIT_RE.split(text) if part)


def _basis_mentions(text: str, catalog: dict[str, Any]) -> list[str]:
    lowered = text or ""
    found: list[str] = []
    seen: set[str] = set()
    display_names = catalog.get("display_name_to_key", {})
    for display_name in sorted(display_names, key=len, reverse=True):
        if _display_name_in_text(display_name, lowered):
            key = display_names[display_name]
            if key not in seen:
                found.append(display_name)
                seen.add(key)
    return found


def _display_name_in_text(display_name: str, text: str) -> bool:
    pattern = re.escape(display_name).replace(r"\-", r"[-\s]?")
    return (
        re.search(
            rf"(?<![A-Za-z0-9]){pattern}(?![A-Za-z0-9])",
            text,
            re.IGNORECASE,
        )
        is not None
    )


def _near_matches(
    normalized: str, catalog: dict[str, Any], limit: int = 6
) -> tuple[str, ...]:
    if not normalized:
        return ()
    aliases = catalog.get("aliases", {})
    scored = []
    for alias, key in aliases.items():
        if normalized in alias or alias in normalized:
            display = catalog["basis_sets"][key]["display_name"]
            scored.append((abs(len(alias) - len(normalized)), display))
    return tuple(name for _, name in sorted(scored)[:limit])


def _family_candidates(
    text: str, catalog: dict[str, Any], limit: int = 8
) -> list[str]:
    lowered = (text or "").lower()
    family = None
    if "karlsruhe" in lowered or "def2" in lowered or "ahlrichs" in lowered:
        family = "ahlrichs"
    elif "dunning" in lowered or "cc-" in lowered:
        family = "dunning"
    elif "pople" in lowered or "31g" in lowered:
        family = "pople"
    if family is None:
        return []
    names = [
        entry["display_name"]
        for entry in catalog.get("basis_sets", {}).values()
        if entry.get("family") == family and entry.get("role") == "orbital"
    ]
    return sorted(names)[:limit]


def _expand_common_spoken_forms(text: str) -> str:
    lowered = (text or "").lower()
    replacements = {
        "six thirty one": "6-31",
        "six three one": "6-31",
        "six dash thirty one": "6-31",
        "six thirty-one": "6-31",
        "double star": "**",
        "star star": "**",
        "single star": "*",
        "tee zeta": "tz",
        "three zeta": "tz",
    }
    for old, new in replacements.items():
        lowered = lowered.replace(old, new)
    return lowered


def _query_intent(query: str) -> dict[str, str | None | set[str]]:
    lowered = _expand_common_spoken_forms(query)
    family = None
    for family_name, terms in _FAMILY_TERMS.items():
        if any(term in lowered for term in terms):
            family = family_name
            break

    quality = {
        label
        for label, terms in _QUALITY_TERMS.items()
        if any(term in lowered for term in terms)
    }
    return {"family": family, "quality": quality}


def _infer_role(query: str, requested: BasisRole) -> BasisRole:
    if requested != "any":
        return requested
    lowered = _expand_common_spoken_forms(query)
    for role, terms in _ROLE_TERMS.items():
        if any(term in lowered for term in terms):
            return role  # type: ignore[return-value]
    return "any"


def _score_entry(
    *,
    key: str,
    entry: dict[str, Any],
    normalized_query: str,
    query_intent: dict[str, str | None | set[str]],
    requested_role: BasisRole,
) -> tuple[int, list[str]]:
    display = entry.get("display_name", "")
    normalized_display = normalize_basis_name(display)
    score, reasons = _name_match_score(normalized_query, normalized_display)

    query_tokens = set(
        _TOKEN_SPLIT_RE.split(_expand_common_spoken_forms(normalized_query))
    )
    display_tokens = set(_TOKEN_SPLIT_RE.split(normalized_display))
    overlap = {token for token in query_tokens & display_tokens if token}
    if overlap:
        score += 12 * len(overlap)
        reasons.append("token_overlap")

    family = query_intent.get("family")
    if family and entry.get("family") == family:
        score += 40
        reasons.append(f"family:{family}")

    quality = query_intent.get("quality") or set()
    if isinstance(quality, set):
        quality_score, quality_reasons = _quality_score(
            normalized_display,
            quality,
        )
        score += quality_score
        reasons.extend(quality_reasons)

    role = entry.get("role")
    if requested_role != "any" and role == requested_role:
        score += 55
        reasons.append(f"role:{requested_role}")
    elif requested_role == "any" and role == "orbital":
        score += 8
        reasons.append("role:orbital")

    if score == 8 and role == "orbital":
        return 0, []
    if key in {"sto-3g", "3-21g"} and "cheap" not in normalized_query:
        score -= 10
    return score, reasons


def _name_match_score(
    normalized_query: str,
    normalized_display: str,
) -> tuple[int, list[str]]:
    score = 0
    reasons: list[str] = []
    if not normalized_query:
        return score, reasons
    if normalized_query == normalized_display:
        score += 200
        reasons.append("exact_name")
    elif normalized_query in normalized_display:
        score += 90
        reasons.append("name_contains_query")
    elif normalized_display in normalized_query:
        score += 80
        reasons.append("query_contains_name")
    if "631" in normalized_query and normalized_display.startswith("631g"):
        score += 70
        reasons.append("pople_shorthand:6-31")
    if normalized_query.startswith("631") and normalized_display == "631g*":
        score += 80
        reasons.append("spoken_name:six_thirty_one_star")
    if "def2tzvp" in normalized_query and "def2tzvp" in normalized_display:
        score += 90
        reasons.append("base_basis:def2-tzvp")
    return score, reasons


def _quality_score(
    normalized_display: str,
    quality: set[str],
) -> tuple[int, list[str]]:
    score = 0
    reasons: list[str] = []
    if "double_zeta" in quality and any(
        tag in normalized_display for tag in ("svp", "dz")
    ):
        score += 26
        reasons.append("quality:double_zeta")
    if "triple_zeta" in quality and any(
        tag in normalized_display for tag in ("tzvp", "tzv", "tz")
    ):
        score += 30
        reasons.append("quality:triple_zeta")
    if "quadruple_zeta" in quality and any(
        tag in normalized_display for tag in ("qzvp", "qzv", "qz")
    ):
        score += 30
        reasons.append("quality:quadruple_zeta")
    if "diffuse" in quality and any(
        tag in normalized_display
        for tag in ("aug", "+", "svpd", "tzvpd", "qzvpd")
    ):
        score += 26
        reasons.append("quality:diffuse")
    if "polarized" in quality and any(
        tag in normalized_display for tag in ("*", "p", "d")
    ):
        score += 18
        reasons.append("quality:polarized")
    return score, reasons


def _candidate_payload(
    entry: dict[str, Any],
    score: int,
    reasons: list[str],
) -> dict[str, Any]:
    payload = {
        "name": entry.get("display_name"),
        "family": entry.get("family"),
        "role": entry.get("role"),
        "score": score,
        "match_reason": reasons[:5],
        "elements_count": len(entry.get("elements", [])),
    }
    auxiliaries = entry.get("auxiliaries") or {}
    if auxiliaries:
        payload["auxiliaries"] = dict(sorted(auxiliaries.items())[:4])
    return payload


def _qualitative_intent_only(
    intent: dict[str, str | None | set[str]],
    normalized_query: str,
) -> bool:
    return bool(intent.get("family") or intent.get("quality")) and not any(
        token in normalized_query
        for token in ("631", "def2", "ccp", "augcc", "lanl", "sto")
    )


def _search_message(
    verdict: BasisIntentVerdict,
    candidates: list[dict[str, Any]],
) -> str:
    if not candidates:
        return "no BSE-backed basis-set candidates matched the query"
    if verdict == "ask_user":
        return "qualitative basis intent matched candidates; ask user or choose conservatively with evidence"
    if verdict == "warn":
        return "multiple close basis-set candidates matched; preserve ambiguity in the answer"
    return "basis-set search returned ranked BSE-backed candidates"
