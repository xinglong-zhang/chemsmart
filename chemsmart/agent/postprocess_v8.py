#!/usr/bin/env python3
"""Deterministic v8 output postprocessor — applied to the model's NORMAL (unconstrained) output before
the runtime consumes it. Uses the chemsmart module index (v8_kind_index.KIND_SETTINGS) to guarantee every
job's `settings` is module-valid, WITHOUT the masking-drop pathology of constrained decoding.

Fixes (the dominant failure classes from the eval):
  * strip per-kind-INVALID settings keys (e.g. `freq` on nci/wbi, `recalc_hess` on orca.opt, the loose
    `calcfc`/`noeigentest`/`ts` keys, a nested `ts:{...}` blob, and any runtime-owned key);
  * reconstruct the TS route into the canonical LIST `additional_opt_options_in_route:["ts","calcfc",...]`
    from the mangled forms the model emits (a comma-string, loose booleans, or a nested `ts` dict);
  * coerce other array-typed keys from string -> list; drop `freq:false` (the default).
It does NOT change kind/file/charge/mult and cannot invent a missing setting or fix a wrong kind — those
are model-quality issues. Importable + self-test."""
from __future__ import annotations
import re

try:
    from . import v8_kind_index as KI
except Exception:  # standalone use outside the package
    try:
        import v8_kind_index as KI
    except Exception:
        KI = None

ROUTE = ("ts", "calcfc", "noeigentest")
_ARRAY_KEYS = {"additional_opt_options_in_route"}


def _split(s):
    return [x for x in re.split(r"[,\s]+", str(s).strip()) if x]


def _route_opts(se, is_ts=False):
    """Recover the TS route list from every mangled form the model emits: a proper list, a comma-string,
    a nested `ts:{...}` dict, a string-valued `ts:"calcfc,noeigentest"` key, or loose `calcfc:true`
    booleans. For a *.ts kind, ensure the leading `ts` token is present."""
    opts = []
    aor = se.get("additional_opt_options_in_route")
    if isinstance(aor, str):
        opts += _split(aor)
    elif isinstance(aor, list):
        opts += [str(x) for x in aor]
    for key in ("ts", "calcfc", "noeigentest"):
        v = se.get(key, None)
        if v is None and key not in se:
            continue
        opts.append(key)                      # the key itself is a route token
        if isinstance(v, str):
            opts += _split(v)
        elif isinstance(v, dict):
            a2 = v.get("additional_opt_options_in_route")
            opts += _split(a2) if isinstance(a2, str) else ([str(x) for x in a2] if isinstance(a2, list) else [])
            opts += [t for t in ROUTE if v.get(t) is True]
    seen, out = set(), []
    for o in opts:
        o = str(o).strip().lower()
        if o and o not in seen:
            seen.add(o)
            out.append(o)
    if is_ts and out and "ts" not in out:    # a TS job's route must start with `ts`
        out = ["ts"] + out
    return out


def postprocess(spec):
    if not isinstance(spec, dict) or spec.get("intent") != "workflow" or KI is None:
        return spec
    for job in spec.get("jobs", []):
        kind = job.get("kind")
        allowed = set(KI.KIND_SETTINGS.get(kind, {}))
        se = dict(job.get("settings", {}) or {})
        new = {}
        # 1) TS route reconstruction (only if the kind accepts it)
        if isinstance(kind, str) and kind.endswith(".ts") and "additional_opt_options_in_route" in allowed:
            r = _route_opts(se, is_ts=True)
            if r:
                new["additional_opt_options_in_route"] = r
        # 2) keep only module-accepted keys; coerce + canonicalize
        for k, v in se.items():
            if k not in allowed or k in new:
                continue
            if k in _ARRAY_KEYS and isinstance(v, str):
                v = _split(v)
            if k == "freq" and v is False:
                continue  # default; drop
            new[k] = v
        if new:
            job["settings"] = new
        elif "settings" in job:
            del job["settings"]
    return spec


def _selftest():
    import copy
    cases = [
        # (name, kind, in_settings, expected_out_settings)
        ("nested ts blob", "gaussian.ts",
         {"ts": {"calcfc": True, "noeigentest": True, "additional_opt_options_in_route": "ts,calcfc,noeigentest"}, "freq": True},
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"], "freq": True}),
        ("loose calcfc/noeig (+ts)", "gaussian.ts", {"calcfc": True, "noeigentest": True, "freq": True},
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"], "freq": True}),
        ("string-valued ts key", "gaussian.ts", {"ts": "calcfc,noeigentest", "freq": True},
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"], "freq": True}),
        ("route as string", "gaussian.ts", {"additional_opt_options_in_route": "ts,calcfc,noeigentest", "freq": True},
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"], "freq": True}),
        ("freq on nci", "gaussian.nci", {"freq": True}, None),
        ("recalc_hess on opt", "orca.opt", {"recalc_hess": True}, None),
        ("valid tddft kept", "gaussian.tddft", {"nstates": 3}, {"nstates": 3}),
        ("runtime-owned stripped", "gaussian.opt", {"functional": "B3LYP", "freq": True}, {"freq": True}),
        ("already-correct ts", "gaussian.ts", {"freq": True, "additional_opt_options_in_route": ["ts", "calcfc"]},
         {"freq": True, "additional_opt_options_in_route": ["ts", "calcfc"]}),
    ]
    ok = 0
    for name, kind, sin, sout in cases:
        spec = {"intent": "workflow", "jobs": [{"id": 1, "kind": kind, "file": "x.xyz", "charge": 0, "mult": 1, "settings": copy.deepcopy(sin)}]}
        postprocess(spec)
        got = spec["jobs"][0].get("settings")
        good = got == sout
        ok += good
        print(f"  {'OK ' if good else 'BAD'} {name:24s} -> {got}")
    print(f"\npostprocessor self-test: {ok}/{len(cases)}")


if __name__ == "__main__":
    _selftest()
