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

# --- v9 hardening (verified +8/45, 0 regressions on the easy/diverse study; deployable, no retrain) ---
# The model HALLUCINATES default-valued / extra keys the clean training data never carried. Strip them at
# the source so the canonical command emerges.
_DEFAULT_VALUED = {"orca.ts": {"tssearch_type": "tssearch"}}   # drop key when its value == the CLI default
_DROP_KEYS = {"gaussian.tddft": ("states", "root", "eqsolv")}  # canonical tddft is nstates-only
_FORCE_PRESENT = {"gaussian.ts": {"freq": True}}               # a Gaussian TS implies freq; add iff omitted


def _split(s):
    return [x for x in re.split(r"[,\s]+", str(s).strip()) if x]


def _fix_scan_type(sd):
    """gaussian.scan: the coordinate-type letter must match the atom count (B=2, A=3, D=4). The model
    sometimes picks the wrong letter (e.g. `B 1 2 3` for a 3-atom angle) -> repair it from the count."""
    if not isinstance(sd, str):
        return sd
    toks = sd.split()
    if len(toks) >= 2 and toks[0] in ("B", "A", "D") and "S" in toks:
        natoms = toks.index("S") - 1
        correct = {2: "B", 3: "A", 4: "D"}.get(natoms)
        if correct and correct != toks[0]:
            toks[0] = correct
            return " ".join(toks)
    return sd


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
        # v9.0) drop hallucinated default-valued / noise keys at the source
        for k, dv in _DEFAULT_VALUED.get(kind, {}).items():
            if se.get(k) == dv:
                se.pop(k, None)
        for k in _DROP_KEYS.get(kind, ()):
            se.pop(k, None)
        if kind == "gaussian.scan" and "scan_definition" in se:
            se["scan_definition"] = _fix_scan_type(se["scan_definition"])
        new = {}
        # 1) TS route reconstruction (only if the kind accepts it)
        if isinstance(kind, str) and kind.endswith(".ts") and "additional_opt_options_in_route" in allowed:
            r = _route_opts(se, is_ts=True)
            if not r and kind == "gaussian.ts":
                r = list(ROUTE)  # v9: inject the canonical route when the model emitted none
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
        # v9.3) ensure canonical-present keys the model omitted (e.g. a Gaussian TS implies freq);
        # only when the key was absent entirely (an explicit False is respected/dropped above)
        for k, v in _FORCE_PRESENT.get(kind, {}).items():
            if k in allowed and k not in new and k not in se:
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
        # --- v9 hardening ---
        ("orca.ts default tssearch stripped", "orca.ts", {"tssearch_type": "tssearch"}, None),
        ("tddft hallucinated keys dropped", "gaussian.tddft",
         {"nstates": 3, "states": "1,2,3", "root": "1", "eqsolv": "HF"}, {"nstates": 3}),
        ("gaussian.ts freq forced when omitted", "gaussian.ts",
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"]},
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"], "freq": True}),
        ("gaussian.ts route forced when omitted", "gaussian.ts", {"freq": True},
         {"freq": True, "additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"]}),
        ("scan wrong type letter fixed", "gaussian.scan", {"scan_definition": "B 1 2 3 S 8 5.0"},
         {"scan_definition": "A 1 2 3 S 8 5.0"}),
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
