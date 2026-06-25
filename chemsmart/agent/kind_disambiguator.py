#!/usr/bin/env python3
"""Deterministic kind-disambiguation backstop for the chemsmart agent  — the lever after the data lever failed TWICE
(v10 templated + v11 diverse contrast both left modred 0/3, wbi 2/3). Operates on (user_query, spec) and
fixes the 3 confusions the 3B prior won't unlearn, using STRONG query keywords + the atoms the model already
extracted. Conservative: only fires on unambiguous keywords; never touches a genuine opt/dias.

Rules:
  R1 wbi: query has Wiberg/bond-index/bond-order AND kind is *.dias|*.wbi -> force *.wbi (drop fragments).
          ALSO spec-only safety: a *.dias with NO fragment_indices is invalid -> treat as *.wbi.
  R2 modred-from-opt: query has modredundant / freeze|constrain|fix|hold the BOND|DISTANCE|ANGLE|DIHEDRAL
          AND kind is *.opt carrying freeze_atoms -> convert to *.modred with modred=[the freeze_atoms].
  R3 modred coords: a *.modred with no coordinates -> extract atom groups from the query.
"""
import re

_WIBERG = re.compile(r"\bwiberg\b|bond[\s-]*index|bond[\s-]*order", re.I)
_INTERNAL = re.compile(r"modredundant"
                       r"|\b(freez\w*|constrain\w*|fix\w*|hold\w*|keep\w*|lock\w*)\b[^.]*\b(bond|distance|angle|dihedral|length)\b"
                       r"|\b(bond|distance|angle|dihedral|length)\b[^.]*\b(fixed|frozen|constrained|held|kept)\b", re.I)
_PAIR = re.compile(r"(\d+)\s*[-–to&,]+\s*(\d+)")


def _atoms_from_freeze(fa):
    if isinstance(fa, str):
        return [int(x) for x in re.split(r"[,\s]+", fa.strip()) if x.lstrip("-").isdigit()]
    if isinstance(fa, list):
        return [int(x) for x in fa]
    return []


def disambiguate(query, spec):
    """Return (spec, changed_bool). spec is mutated in place."""
    if not isinstance(spec, dict) or spec.get("intent") != "workflow":
        return spec, False
    q = query or ""
    changed = False
    for job in spec.get("jobs", []):
        k = job.get("kind") or ""
        s = job.get("settings") or {}
        prog = k.split(".")[0] if "." in k else "gaussian"
        # R1 — Wiberg / bond-index/order, or an invalid fragment-less dias
        if k.endswith(".dias") and (_WIBERG.search(q) or "fragment_indices" not in s):
            job["kind"] = f"{prog}.wbi"; job.pop("settings", None); changed = True; continue
        if k.endswith(".wbi") and _WIBERG.search(q):
            job.pop("settings", None)  # canonical wbi has no settings
        # R2 — internal-coordinate constraint mislabeled as opt+freeze_atoms
        if k.endswith(".opt") and s.get("freeze_atoms") and _INTERNAL.search(q):
            atoms = _atoms_from_freeze(s["freeze_atoms"])
            if len(atoms) >= 2:
                job["kind"] = f"{prog}.modred"
                job["settings"] = {"modred": [atoms]}
                changed = True; continue
        # R3 — modred with no coordinates: recover atom groups from the query
        if k.endswith(".modred") and not s.get("modred"):
            pairs = [[int(a), int(b)] for a, b in _PAIR.findall(q)]
            if pairs:
                s = dict(s); s["modred"] = pairs; job["settings"] = s; changed = True
    return spec, changed
