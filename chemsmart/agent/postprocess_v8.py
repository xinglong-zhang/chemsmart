#!/usr/bin/env python3
"""Deterministic v8 output postprocessor — applied to the model's NORMAL (unconstrained) output before
the runtime consumes it. Uses the chemsmart module index (v8_kind_index.KIND_SETTINGS) to guarantee every
job's `settings` is module-valid, WITHOUT the masking-drop pathology of constrained decoding.

Fixes (the dominant failure classes from the eval):
  * strip per-kind-INVALID settings keys (e.g. `freq` on nci/wbi, `recalc_hess` on orca.opt, the loose
    `calcfc`/`noeigentest`/`ts` keys, a nested `ts:{...}` blob, and any runtime-owned key);
  * reconstruct the TS route into the canonical LIST `additional_opt_options_in_route:["ts","calcfc",...]`
    from the mangled forms the model emits (a comma-string, loose booleans, or a nested `ts` dict);
  * coerce other array-typed keys from string -> list; drop `freq:false` (the default);
  * canonicalize the known synthetic `*.opt+freq` / `*.opt.freq` kinds into `*.opt` with `freq:true`.
It does NOT change file/charge/mult and cannot invent missing structural settings — those are model-quality
issues. Importable + self-test."""
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

# --- v9 hardening (verified against the real chemsmart CLI; diversity-preserving) ---
# Rules are limited to what the repo scan confirmed is safe. We do NOT strip `states/root/eqsolv` from
# gaussian.tddft: the repo (cli/gaussian/tddft.py `def td(ctx, states, root, nstates, eqsolv, ...)`) shows
# these are REAL, valid TD-DFT options — stripping them would destroy legitimate requests (e.g. root=2).
# The model over-emitting derived values is a training-side issue, not a postprocess one.
#
# orca.ts: `--tssearch-type` is real (cli/orca/ts.py) with default "optts" and values {optts, scants}. The
# model sometimes emits the bogus value "tssearch" or the redundant default "optts"; drop either (-> the
# runtime default), but KEEP a meaningful non-default like "scants".
_TSSEARCH_DROP = {"tssearch", "optts"}
# B1 verified: bare `gaussian ts` already emits opt=(ts,calcfc,noeigentest) (calcfc = the Hessian); a TS
# does NOT auto-imply a separate freq job. So do not force freq — the model emits freq only if requested.
_FORCE_PRESENT = {}
_KIND_CANON = {
    "gaussian.opt+freq": "gaussian.opt",
    "gaussian.opt.freq": "gaussian.opt",
    "gaussian.opt_freq": "gaussian.opt",
    "orca.opt+freq": "orca.opt",
    "orca.opt.freq": "orca.opt",
    "orca.opt_freq": "orca.opt",
}


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
        raw_kind = str(kind).strip().lower().replace(" ", "")
        canon_kind = _KIND_CANON.get(raw_kind)
        se = dict(job.get("settings", {}) or {})
        if canon_kind:
            job["kind"] = canon_kind
            kind = canon_kind
            se["freq"] = True
        allowed = set(KI.KIND_SETTINGS.get(kind, {}))
        # v9.0) drop the redundant-default / bogus orca.ts search type (keeps meaningful values e.g. scants)
        if kind == "orca.ts" and str(se.get("tssearch_type", "")).lower() in _TSSEARCH_DROP:
            se.pop("tssearch_type", None)
        if kind == "gaussian.scan" and "scan_definition" in se:
            se["scan_definition"] = _fix_scan_type(se["scan_definition"])
        new = {}
        # 1) TS route: chemsmart AUTO-derives opt=(ts,calcfc,noeigentest) from a *.ts job (verified B1).
        # Emitting the canonical route duplicates it AND renders the list as a literal -> BROKEN. So keep
        # ONLY genuine EXTRA opts beyond the canonical triple; omit otherwise (let the runtime derive it).
        if isinstance(kind, str) and kind.endswith(".ts") and "additional_opt_options_in_route" in allowed:
            extras = [t for t in _route_opts(se, is_ts=False) if t not in ROUTE]
            if extras:
                new["additional_opt_options_in_route"] = extras
        # 2) keep only module-accepted keys; coerce + canonicalize
        _ts = isinstance(kind, str) and kind.endswith(".ts")
        for k, v in se.items():
            if k not in allowed or k in new:
                continue
            if k == "additional_opt_options_in_route" and _ts:
                continue  # handled in step 1 (extras only; canonical route is auto-derived)
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
        # --- B2: the canonical ts route is auto-derived -> DROP it (any mangled form); keep freq if present ---
        ("ts nested blob -> route dropped", "gaussian.ts",
         {"ts": {"calcfc": True, "noeigentest": True, "additional_opt_options_in_route": "ts,calcfc,noeigentest"}, "freq": True},
         {"freq": True}),
        ("ts loose calcfc/noeig -> dropped", "gaussian.ts", {"calcfc": True, "noeigentest": True, "freq": True},
         {"freq": True}),
        ("ts string-valued key -> dropped", "gaussian.ts", {"ts": "calcfc,noeigentest", "freq": True},
         {"freq": True}),
        ("ts route-as-string -> dropped", "gaussian.ts", {"additional_opt_options_in_route": "ts,calcfc,noeigentest", "freq": True},
         {"freq": True}),
        ("ts canonical (no freq) -> empty", "gaussian.ts",
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest"]}, None),
        ("ts EXTRA opt kept", "gaussian.ts",
         {"additional_opt_options_in_route": ["ts", "calcfc", "noeigentest", "maxstep=15"]},
         {"additional_opt_options_in_route": ["maxstep=15"]}),
        ("freq on nci", "gaussian.nci", {"freq": True}, None),
        ("recalc_hess on opt", "orca.opt", {"recalc_hess": True}, None),
        ("valid tddft kept", "gaussian.tddft", {"nstates": 3}, {"nstates": 3}),
        ("runtime-owned stripped", "gaussian.opt", {"functional": "B3LYP", "freq": True}, {"freq": True}),
        ("orca.opt freq KEPT (B2)", "orca.opt", {"freq": True}, {"freq": True}),
        ("gaussian.opt+freq kind canonicalized", "gaussian.opt+freq", {}, {"freq": True}),
        ("orca.opt.freq kind canonicalized", "orca.opt.freq", {"additional_route_parameters": "tightscf"},
         {"additional_route_parameters": "tightscf", "freq": True}),
        # --- v9 hardening (CLI-verified, diversity-preserving) ---
        ("orca.ts bogus tssearch dropped", "orca.ts", {"tssearch_type": "tssearch"}, None),
        ("orca.ts default optts dropped", "orca.ts", {"tssearch_type": "optts"}, None),
        ("orca.ts scants KEPT (diversity)", "orca.ts", {"tssearch_type": "scants"}, {"tssearch_type": "scants"}),
        ("tddft real options KEPT (diversity)", "gaussian.tddft",
         {"nstates": 3, "root": 2}, {"nstates": 3, "root": 2}),
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
