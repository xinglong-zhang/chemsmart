#!/usr/bin/env python3
"""v8 MODULE INDEX: per-kind model-emittable `settings` keys, derived from the chemsmart repo's actual
CLI subcommand options (chemsmart/cli/{gaussian,orca}/<kind>.py) + the per-kind settings subclasses,
MINUS the runtime-owned keys. This is authoritative from the code — independent of any eval's expected
answer. Used to (a) judge whether a model answer is a VALID chemsmart command for its kind, and
(b) build the per-kind JSON-schema grammar for constrained decoding.

Each key maps to a JSON-schema type tag (so the grammar can also fix structure, e.g. force
additional_opt_options_in_route to be an ARRAY, not the string the model sometimes emits)."""

# type tags -> JSON schema fragment
_T = {
    "bool": {"type": "boolean"},
    "int": {"type": "integer"},
    "num": {"type": "number"},
    "str": {"type": "string"},
    "dir": {"type": "string", "enum": ["forward", "backward", "both", "reverse"]},
    "strs": {"type": "array", "items": {"type": "string"}},
    "ints": {"type": "array", "items": {"type": "integer"}},
    "pairs": {"type": "array", "items": {"type": "array", "items": {"type": "integer"}}},
    "any": {},
}

# kind -> {settings_key: type_tag}. Empty dict = NO settings allowed for that kind.
KIND_SETTINGS = {
    "gaussian.opt":   {"freq": "bool", "freeze_atoms": "ints", "additional_route_parameters": "str"},
    "gaussian.sp":    {},
    "gaussian.freq":  {},
    "gaussian.ts":    {"freq": "bool", "additional_opt_options_in_route": "strs", "freeze_atoms": "ints",
                       "additional_route_parameters": "str"},
    "gaussian.irc":   {"direction": "dir", "flat_irc": "bool", "additional_route_parameters": "str"},
    "gaussian.scan":  {"scan_definition": "str", "additional_route_parameters": "str"},
    "gaussian.modred": {"modred": "pairs", "additional_route_parameters": "str"},
    "gaussian.nci":   {},
    "gaussian.resp":  {},
    "gaussian.tddft": {"nstates": "int", "states": "any", "root": "any", "eqsolv": "any",
                       "additional_route_parameters": "str"},
    "gaussian.dias":  {"fragment_indices": "ints"},  # fragment 1 ONLY (runtime derives fragment 2 complement)
    "gaussian.crest": {"num_confs_to_run": "int", "grouping_strategy": "str", "num_groups": "int"},
    "gaussian.traj":  {"num_structures_to_run": "int", "proportion_structures_to_use": "num",
                       "grouping_strategy": "str"},
    "gaussian.wbi":   {},
    "gaussian.qmmm":  {"high_level_atoms": "ints", "low_level_atoms": "ints"},
    "gaussian.qrc":   {},
    "orca.sp":        {},
    "orca.opt":       {"freq": "bool", "freeze_atoms": "ints", "additional_route_parameters": "str"},
    "orca.ts":        {"freq": "bool", "recalc_hess": "bool", "trust_radius": "num", "tssearch_type": "str",
                       "additional_route_parameters": "str"},
    "orca.freq":      {},
    "orca.irc":       {"direction": "dir", "inithess": "str", "additional_route_parameters": "str"},
    "orca.scan":      {"scan_definition": "str", "additional_route_parameters": "str"},
    "orca.modred":    {"modred": "pairs", "additional_route_parameters": "str"},
    "orca.neb":       {"nimages": "int", "joboption": "str"},
    "orca.qmmm":      {"high_level_atoms": "ints"},
    "orca.qrc":       {},
}
# `additional_route_parameters` is a GROUP-LEVEL gaussian/orca CLI option (cli/gaussian/gaussian.py:578,
# cli/orca/orca.py:444) -> valid on EVERY subcommand. Allow it on all kinds so the model can pass a real
# route keyword (e.g. scf=tight) on sp/freq/etc. and the postprocessor does not silently drop it.
for _kind in KIND_SETTINGS:
    KIND_SETTINGS[_kind].setdefault("additional_route_parameters", "str")

GEOM_KINDS = {"opt", "ts", "irc", "scan", "modred", "crest", "neb"}


def allowed(kind):
    return set(KIND_SETTINGS.get(kind, {}))


def settings_keys_valid(kind, settings):
    """True iff every emitted settings key is module-accepted for this kind."""
    return all(k in KIND_SETTINGS.get(kind, {}) for k in (settings or {}))


def _job_schema(kind):
    keys = KIND_SETTINGS[kind]
    settings = {"type": "object", "additionalProperties": False,
                "properties": {k: _T[t] for k, t in keys.items()}}
    props = {
        "id": {"type": "integer"},
        "kind": {"const": kind},
        "charge": {"type": "integer"},
        "mult": {"type": "integer"},
        "file": {"type": "string"},
        "geom_from": {"type": "integer"},
        "label": {"type": "string"},
        "execution": {"type": "string", "enum": ["run_local", "submit"]},
        "server": {"type": "string"},
    }
    if keys:
        props["settings"] = settings
    if kind == "orca.neb":
        props["product_file"] = {"type": "string"}
    # the molecule input is mandatory — REQUIRE it, else the constrained model drops `file` entirely.
    # (single-job study uses `file`; chains with geom_from need a generalized schema later.)
    req = ["id", "kind", "charge", "mult", "file"]
    return {"type": "object", "additionalProperties": False, "required": req, "properties": props}


def build_spec_schema():
    """The full v8 output JSON schema: per-kind discriminated jobs (settings constrained per kind) for
    workflow, or {intent, message} for non-workflow."""
    workflow = {
        "type": "object", "additionalProperties": False,
        "required": ["intent", "jobs"],
        "properties": {
            "intent": {"const": "workflow"},
            "jobs": {"type": "array", "minItems": 1,
                     "items": {"oneOf": [_job_schema(k) for k in KIND_SETTINGS]}},
        },
    }
    nonworkflow = {
        "type": "object", "additionalProperties": False,
        "required": ["intent", "message"],
        "properties": {"intent": {"type": "string", "enum": ["advisory", "decline", "chitchat"]},
                       "message": {"type": "string"}},
    }
    return {"oneOf": [workflow, nonworkflow]}


if __name__ == "__main__":
    import json
    print(f"kinds: {len(KIND_SETTINGS)}")
    print("schema size (chars):", len(json.dumps(build_spec_schema())))
    # sanity: schema is valid JSON-schema-ish
    s = build_spec_schema()
    assert "oneOf" in s
    print("OK")
