"""Deterministic cleanup for compact chemsmart agent SPEC output."""

from __future__ import annotations

import re
from typing import Any

from chemsmart.agent import v8_kind_index as KI

ROUTE = ("ts", "calcfc", "noeigentest")
_ARRAY_KEYS = {"additional_opt_options_in_route"}
_TSSEARCH_DROP = {"tssearch", "optts"}
_TDDFT_STATES = {"singlets", "triplets", "50-50"}
_KIND_CANON = {
    "gaussian.opt+freq": "gaussian.opt",
    "gaussian.opt.freq": "gaussian.opt",
    "gaussian.opt_freq": "gaussian.opt",
    "orca.opt+freq": "orca.opt",
    "orca.opt.freq": "orca.opt",
    "orca.opt_freq": "orca.opt",
}
_ALLOWED_TS_EXTRA = re.compile(
    r"^(calcall|cartesian|tight|"
    r"(maxstep|maxcycles|maxcycle|recalcfc|maxmicroiterations)=\d+)$",
    re.IGNORECASE,
)


def _split(value: Any) -> list[str]:
    return [
        item
        for item in re.split(r"[,\s/]+", str(value).strip())
        if item
    ]


def _fix_scan_type(scan_definition: Any) -> Any:
    if not isinstance(scan_definition, str):
        return scan_definition
    tokens = scan_definition.split()
    if len(tokens) >= 2 and tokens[0] in {"B", "A", "D"} and "S" in tokens:
        atom_count = tokens.index("S") - 1
        correct = {2: "B", 3: "A", 4: "D"}.get(atom_count)
        if correct and correct != tokens[0]:
            tokens[0] = correct
            return " ".join(tokens)
    return scan_definition


def _route_opts(settings: dict[str, Any], *, is_ts: bool = False) -> list[str]:
    opts: list[str] = []
    route = settings.get("additional_opt_options_in_route")
    if isinstance(route, str):
        opts += _split(route)
    elif isinstance(route, list):
        for item in route:
            opts += _split(item)
    for key in ROUTE:
        if key not in settings:
            continue
        value = settings.get(key)
        opts.append(key)
        if isinstance(value, str):
            opts += _split(value)
        elif isinstance(value, dict):
            nested = value.get("additional_opt_options_in_route")
            if isinstance(nested, str):
                opts += _split(nested)
            elif isinstance(nested, list):
                opts += [str(item) for item in nested]
            opts += [token for token in ROUTE if value.get(token) is True]
    seen: set[str] = set()
    out: list[str] = []
    for option in opts:
        option = str(option).strip().lower()
        if option and option not in seen:
            seen.add(option)
            out.append(option)
    if is_ts and out and "ts" not in out:
        out.insert(0, "ts")
    return out


def postprocess(spec: dict[str, Any]) -> dict[str, Any]:
    """Normalize compact SPEC before adapter rendering.

    This strips runtime-owned route tokens that chemsmart derives itself,
    keeps genuine per-kind structural settings, and fixes the known
    opt+freq and scan-definition drift classes observed in v11-v13.1 evals.
    """
    if not isinstance(spec, dict) or spec.get("intent") != "workflow":
        return spec
    jobs = spec.get("jobs")
    if not isinstance(jobs, list):
        return spec

    for job in jobs:
        if not isinstance(job, dict):
            continue
        raw_kind = str(job.get("kind", "")).strip().lower().replace(" ", "")
        kind = _KIND_CANON.get(raw_kind, raw_kind)
        if kind:
            job["kind"] = kind
        settings = dict(job.get("settings", {}) or {})
        if raw_kind in _KIND_CANON:
            settings["freq"] = True

        allowed = set(KI.KIND_SETTINGS.get(kind, {}))
        if kind == "orca.ts" and str(
            settings.get("tssearch_type", "")
        ).lower() in _TSSEARCH_DROP:
            settings.pop("tssearch_type", None)
        if kind == "gaussian.scan" and "scan_definition" in settings:
            settings["scan_definition"] = _fix_scan_type(
                settings["scan_definition"]
            )
        if kind == "gaussian.tddft" and "states" in settings:
            states = str(settings["states"]).strip().lower()
            if states not in _TDDFT_STATES:
                settings.pop("states", None)

        new: dict[str, Any] = {}
        is_ts = isinstance(kind, str) and kind.endswith(".ts")
        if is_ts and "additional_opt_options_in_route" in allowed:
            extras = [
                token
                for token in _route_opts(settings, is_ts=False)
                if token not in ROUTE and _ALLOWED_TS_EXTRA.match(token)
            ]
            if extras:
                new["additional_opt_options_in_route"] = extras

        for key, value in settings.items():
            if key not in allowed or key in new:
                continue
            if key == "additional_opt_options_in_route" and is_ts:
                continue
            if key == "additional_route_parameters" and not isinstance(
                value, str
            ):
                continue
            if key in _ARRAY_KEYS and isinstance(value, str):
                value = _split(value)
            if key == "freq" and value is False:
                continue
            new[key] = value

        if new:
            job["settings"] = new
        else:
            job.pop("settings", None)
    return spec


__all__ = ["postprocess"]
