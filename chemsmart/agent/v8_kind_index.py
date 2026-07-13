"""Compact-SPEC kind index for the local chemsmart agent.

The model emits only structural, per-kind settings. Project defaults such as
functional, basis, solvent, and scheduler policy remain runtime-owned.
"""

from __future__ import annotations

_T = {
    "bool": {"type": "boolean"},
    "int": {"type": "integer"},
    "num": {"type": "number"},
    "str": {"type": "string"},
    "dir": {
        "type": "string",
        "enum": ["forward", "backward", "both", "reverse"],
    },
    "qmmm_parent": {
        "type": "string",
        "enum": ["opt", "sp", "ts", "scan", "modred", "qrc", "neb"],
    },
    "strs": {"type": "array", "items": {"type": "string"}},
    "ints": {"type": "array", "items": {"type": "integer"}},
    "pairs": {
        "type": "array",
        "items": {"type": "array", "items": {"type": "integer"}},
    },
    "any": {},
}

KIND_SETTINGS = {
    "gaussian.opt": {
        "freq": "bool",
        "freeze_atoms": "ints",
        "additional_route_parameters": "str",
    },
    "gaussian.sp": {},
    "gaussian.freq": {},
    "gaussian.ts": {
        "freq": "bool",
        "additional_opt_options_in_route": "strs",
        "freeze_atoms": "ints",
        "additional_route_parameters": "str",
    },
    "gaussian.irc": {
        "direction": "dir",
        "flat_irc": "bool",
        "additional_route_parameters": "str",
    },
    "gaussian.scan": {
        "scan_definition": "str",
        "additional_route_parameters": "str",
    },
    "gaussian.modred": {
        "modred": "pairs",
        "additional_route_parameters": "str",
    },
    "gaussian.nci": {},
    "gaussian.resp": {},
    "gaussian.tddft": {
        "nstates": "int",
        "states": "any",
        "root": "any",
        "eqsolv": "any",
        "additional_route_parameters": "str",
    },
    # Fragment 1 only. The runtime derives fragment 2 as the complement.
    "gaussian.dias": {"fragment_indices": "ints"},
    "gaussian.crest": {
        "num_confs_to_run": "int",
        "grouping_strategy": "str",
        "num_groups": "int",
    },
    "gaussian.traj": {
        "num_structures_to_run": "int",
        "proportion_structures_to_use": "num",
        "grouping_strategy": "str",
    },
    "gaussian.wbi": {},
    "gaussian.qmmm": {
        "parent_job": "qmmm_parent",
        "high_level_atoms": "ints",
        "medium_level_atoms": "ints",
        "low_level_atoms": "ints",
        "charge_total": "int",
        "mult_total": "int",
        "charge_intermediate": "int",
        "mult_intermediate": "int",
        "charge_high": "int",
        "mult_high": "int",
    },
    "gaussian.qrc": {},
    "orca.sp": {},
    "orca.opt": {
        "freq": "bool",
        "freeze_atoms": "ints",
        "additional_route_parameters": "str",
    },
    "orca.ts": {
        "freq": "bool",
        "recalc_hess": "bool",
        "trust_radius": "num",
        "tssearch_type": "str",
        "additional_route_parameters": "str",
    },
    "orca.freq": {},
    "orca.irc": {
        "direction": "dir",
        "inithess": "str",
        "additional_route_parameters": "str",
    },
    "orca.scan": {
        "scan_definition": "str",
        "additional_route_parameters": "str",
    },
    "orca.modred": {
        "modred": "pairs",
        "additional_route_parameters": "str",
    },
    "orca.neb": {"nimages": "int", "joboption": "str"},
    "orca.qmmm": {
        "parent_job": "qmmm_parent",
        "high_level_atoms": "ints",
        "intermediate_level_atoms": "ints",
        "charge_total": "int",
        "mult_total": "int",
        "charge_intermediate": "int",
        "mult_intermediate": "int",
        "charge_high": "int",
        "mult_high": "int",
        "jobtype": "str",
        "low_level_method": "str",
    },
    "orca.qrc": {},
}

# Group-level route parameters are accepted for both Gaussian and ORCA command
# groups, so keep them for every kind instead of silently dropping them.
for _kind in KIND_SETTINGS:
    KIND_SETTINGS[_kind].setdefault("additional_route_parameters", "str")

GEOM_KINDS = {"opt", "ts", "irc", "scan", "modred", "crest", "neb"}


def allowed(kind: str) -> set[str]:
    return set(KIND_SETTINGS.get(kind, {}))


def settings_keys_valid(kind: str, settings: dict | None) -> bool:
    return all(k in KIND_SETTINGS.get(kind, {}) for k in (settings or {}))


def _job_schema(kind: str) -> dict:
    keys = KIND_SETTINGS[kind]
    settings = {
        "type": "object",
        "additionalProperties": False,
        "properties": {k: _T[t] for k, t in keys.items()},
    }
    props = {
        "id": {"type": "integer"},
        "kind": {"const": kind},
        "charge": {"type": "integer"},
        "mult": {"type": "integer"},
        "file": {"type": "string"},
        "geom_from": {"type": "integer"},
        "record_index": {"type": "integer"},
        "record_id": {"type": "string"},
        "structure_index": {"type": "string"},
        "structure_id": {"type": "string"},
        "molecule_id": {"type": "string"},
        "label": {"type": "string"},
        "execution": {"type": "string", "enum": ["run_local", "submit"]},
        "server": {"type": "string"},
    }
    if keys:
        props["settings"] = settings
    if kind == "orca.neb":
        props["product_file"] = {"type": "string"}
    return {
        "type": "object",
        "additionalProperties": False,
        "required": ["id", "kind", "charge", "mult", "file"],
        "properties": props,
    }


def build_spec_schema() -> dict:
    workflow = {
        "type": "object",
        "additionalProperties": False,
        "required": ["intent", "jobs"],
        "properties": {
            "intent": {"const": "workflow"},
            "jobs": {
                "type": "array",
                "minItems": 1,
                "items": {"oneOf": [_job_schema(k) for k in KIND_SETTINGS]},
            },
        },
    }
    nonworkflow = {
        "type": "object",
        "additionalProperties": False,
        "required": ["intent", "message"],
        "properties": {
            "intent": {
                "type": "string",
                "enum": ["advisory", "decline", "chitchat"],
            },
            "message": {"type": "string"},
        },
    }
    return {"oneOf": [workflow, nonworkflow]}
