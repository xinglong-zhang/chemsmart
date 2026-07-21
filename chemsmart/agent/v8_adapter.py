"""Adapter from compact chemsmart SPEC to real CLI commands."""

from __future__ import annotations

import json
import shlex
from typing import Any

from chemsmart.agent.postprocess_v8 import postprocess

# Map SPEC kind suffixes to real CLI subcommands. Gaussian/ORCA single-point is
# `sp` in the real CLI, not `singlepoint`; `*.freq` has no literal subcommand
# and is routed through `opt` with a route-level freq keyword.
_SUBCMD = {"tddft": "td", "freq": "opt"}
_FLAG = {
    "additional_opt_options_in_route": "--additional-opt-options",
    "additional_route_parameters": "--additional-route-parameters",
    "direction": "--direction",
    "flat_irc": "--flat-irc",
    "inithess": "--inithess",
    "maxiter": "--maxiter",
    "hessmode": "--hessmode",
    "nstates": "--nstates",
    "states": "--states",
    "root": "--root",
    "eqsolv": "--eqsolv",
    "fragment_indices": "--fragment-indices",
    "num_confs_to_run": "--num-confs-to-run",
    "grouping_strategy": "--grouping-strategy",
    "num_groups": "--num-groups",
    "num_structures_to_run": "--num-structures-to-run",
    "proportion_structures_to_use": "--proportion-structures-to-use",
    "high_level_atoms": "--high-level-atoms",
    "medium_level_atoms": "--medium-level-atoms",
    "low_level_atoms": "--low-level-atoms",
    "intermediate_level_atoms": "--intermediate-level-atoms",
    "charge_total": "--charge-total",
    "mult_total": "--mult-total",
    "charge_intermediate": "--charge-intermediate",
    "mult_intermediate": "--mult-intermediate",
    "charge_high": "--charge-high",
    "mult_high": "--mult-high",
    "jobtype": "--jobtype",
    "low_level_method": "--low-level-method",
    "nimages": "--nimages",
    "joboption": "--joboption",
    "freeze_atoms": "--freeze-atoms",
    "recalc_hess": "--recalc-hess",
    "trust_radius": "--trust-radius",
    "tssearch_type": "--tssearch-type",
    "scan_definition": "--coordinates",
    "modred": "--coordinates",
}
_RANGE_FLAGS = {
    "--fragment-indices",
    "--high-level-atoms",
    "--medium-level-atoms",
    "--low-level-atoms",
    "--freeze-atoms",
}
_RANGE_SETTINGS = {
    "fragment_indices",
    "high_level_atoms",
    "medium_level_atoms",
    "low_level_atoms",
    "intermediate_level_atoms",
    "freeze_atoms",
}
_DB_SELECTOR_FLAG = {
    "record_index": "--record-index",
    "record_id": "--record-id",
    "structure_index": "--structure-index",
    "structure_id": "--structure-id",
    "molecule_id": "--molecule-id",
}
_SUBCOMMAND_SETTINGS = {
    "fragment_indices",
    "high_level_atoms",
    "medium_level_atoms",
    "low_level_atoms",
    "nimages",
    "joboption",
    "freeze_atoms",
    "recalc_hess",
    "trust_radius",
    "tssearch_type",
    "scan_definition",
    "modred",
    "nstates",
    "states",
    "root",
    "eqsolv",
    "intermediate_level_atoms",
    "charge_total",
    "mult_total",
    "charge_intermediate",
    "mult_intermediate",
    "charge_high",
    "mult_high",
    "jobtype",
    "low_level_method",
}
_QMMM_PARENTS = {
    "gaussian": {"opt", "sp", "ts", "scan", "modred", "qrc"},
    "orca": {"opt", "sp", "ts", "scan", "modred", "qrc", "neb"},
}


def _subcommand(kind: str) -> str:
    suffix = kind.split(".", 1)[1]
    return _SUBCMD.get(suffix, suffix)


def _range_value(value: Any, *, first_fragment_only: bool = False) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        if first_fragment_only and value and isinstance(value[0], list):
            value = value[0]
        flattened: list[Any] = []
        for item in value:
            if isinstance(item, list):
                flattened.extend(item)
            else:
                flattened.append(item)
        return ",".join(str(item) for item in flattened)
    return str(value)


def _pair_value(value: Any) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        groups: list[str] = []
        for item in value:
            if isinstance(item, list):
                groups.append(",".join(str(atom) for atom in item))
            else:
                groups.append(str(item))
        return " ".join(groups)
    return str(value)


def _literal_list_value(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, separators=(",", ":"))


def _parse_scan_definition(scan_definition: str) -> dict[str, list[Any]]:
    lines = [
        line.strip()
        for chunk in str(scan_definition).splitlines()
        for line in chunk.split(";")
        if line.strip()
    ]
    if not lines:
        raise ValueError("scan_definition cannot be empty")

    atom_count_by_coordinate = {"B": 2, "A": 3, "D": 4}
    coords: list[list[int]] = []
    num_steps: list[int] = []
    step_size: list[float] = []
    constrained_coordinates: list[list[int]] = []

    for line in lines:
        tokens = line.split()
        coordinate_type = tokens[0].upper()
        expected_atoms = atom_count_by_coordinate.get(coordinate_type)
        if expected_atoms is None:
            raise ValueError(
                f"scan_definition {line!r} must start with B, A, or D"
            )
        if len(tokens) < expected_atoms + 2:
            raise ValueError(f"scan_definition {line!r} is incomplete")
        try:
            coordinate = [
                int(token) for token in tokens[1 : 1 + expected_atoms]
            ]
        except ValueError as exc:
            raise ValueError(
                f"scan_definition {line!r} has non-integer atom indices"
            ) from exc

        directive = tokens[1 + expected_atoms].upper()
        if directive == "F":
            constrained_coordinates.append(coordinate)
            continue
        if directive != "S" or len(tokens) != expected_atoms + 4:
            raise ValueError(
                "scan_definition entries must look like "
                "'D 1 2 3 4 S 10 36.0' or 'B 1 2 F'"
            )
        try:
            steps = int(tokens[2 + expected_atoms])
            size = float(tokens[3 + expected_atoms])
        except ValueError as exc:
            raise ValueError(
                f"scan_definition {line!r} has invalid steps or step size"
            ) from exc
        coords.append(coordinate)
        num_steps.append(steps)
        step_size.append(size)

    if not coords:
        raise ValueError(
            "scan_definition must include at least one S scan entry"
        )
    parsed: dict[str, list[Any]] = {
        "coords": coords,
        "num_steps": num_steps,
        "step_size": step_size,
    }
    if constrained_coordinates:
        parsed["constrained_coordinates"] = constrained_coordinates
    return parsed


def _fmt_scan_sequence(values: list[Any]) -> str:
    if len(values) == 1:
        return shlex.quote(str(values[0]))
    return shlex.quote(json.dumps(values, separators=(",", ":")))


def _gaussian_scan_flags(scan_definition: Any) -> list[str]:
    parsed = _parse_scan_definition(str(scan_definition))
    flags = [
        "--coordinates",
        shlex.quote(json.dumps(parsed["coords"], separators=(",", ":"))),
        "--num-steps",
        _fmt_scan_sequence(parsed["num_steps"]),
        "--step-size",
        _fmt_scan_sequence(parsed["step_size"]),
    ]
    constrained = parsed.get("constrained_coordinates")
    if constrained:
        flags += [
            "--constrained-coordinates",
            shlex.quote(json.dumps(constrained, separators=(",", ":"))),
        ]
    return flags


def _fmt_setting(key: str, value: Any) -> str:
    if key == "eqsolv" and isinstance(value, bool):
        return "eqsolv" if value else "noneqsolv"
    if isinstance(value, bool):
        return "true" if value else "false"
    if key == "fragment_indices":
        return shlex.quote(_range_value(value, first_fragment_only=True))
    if key in _RANGE_SETTINGS:
        return shlex.quote(_range_value(value))
    if key == "modred":
        return shlex.quote(_literal_list_value(value))
    if key == "scan_definition":
        return shlex.quote(_pair_value(value))
    if isinstance(value, list):
        return shlex.quote(",".join(str(item) for item in value))
    return shlex.quote(str(value))


def _pop_qmmm_parent(
    kind: str, program: str, settings: dict[str, Any]
) -> str | None:
    if not kind.endswith(".qmmm"):
        return None
    qmmm_parent = str(settings.pop("parent_job", "")).strip().lower()
    if qmmm_parent not in _QMMM_PARENTS[program]:
        allowed = ", ".join(sorted(_QMMM_PARENTS[program]))
        raise ValueError(
            f"{kind} requires settings.parent_job; allowed: {allowed}"
        )
    return qmmm_parent


def _settings_flags(
    kind: str, settings: dict[str, Any]
) -> tuple[list[str], list[str]]:
    freq_true = settings.pop("freq", None) is True or kind.endswith(".freq")
    if freq_true:
        route = settings.get("additional_route_parameters")
        settings["additional_route_parameters"] = (
            f"{route} freq" if isinstance(route, str) and route else "freq"
        )
    group_flags: list[str] = []
    subcommand_flags: list[str] = []
    for key, value in settings.items():
        if kind == "gaussian.scan" and key == "scan_definition":
            subcommand_flags += _gaussian_scan_flags(value)
            continue
        flag = _FLAG.get(key)
        if flag:
            target = (
                subcommand_flags
                if key in _SUBCOMMAND_SETTINGS
                else group_flags
            )
            target += [flag, _fmt_setting(key, value)]
    return group_flags, subcommand_flags


def _db_selector_flags(job: dict[str, Any]) -> list[str]:
    flags: list[str] = []
    for key, flag in _DB_SELECTOR_FLAG.items():
        value = job.get(key)
        if value is not None:
            flags += [flag, shlex.quote(str(value))]
    return flags


def _job_command(
    job: dict[str, Any],
    geom_of: dict[Any, str],
    default_project: str | None = None,
) -> str:
    kind = str(job["kind"])
    program = kind.split(".", 1)[0]
    verb = "sub" if job.get("execution") == "submit" else "run"
    parts = ["chemsmart", verb]
    if job.get("server"):
        parts += ["-s", shlex.quote(str(job["server"]))]
    parts.append(program)
    project = job.get("project") or default_project
    if project:
        parts += ["-p", shlex.quote(str(project))]

    settings = dict(job.get("settings", {}) or {})
    qmmm_parent = _pop_qmmm_parent(kind, program, settings)
    group_flags, subcommand_flags = _settings_flags(kind, settings)
    group_flags += _db_selector_flags(job)

    parts += group_flags
    source = job.get("file") or geom_of.get(job.get("geom_from"))
    if not source:
        source = "<upstream-geometry>"
    parts += ["-f", shlex.quote(str(source))]
    if kind == "orca.neb" and job.get("product_file"):
        subcommand_flags += ["-e", shlex.quote(str(job["product_file"]))]
    parts += ["-c", str(job.get("charge", 0)), "-m", str(job.get("mult", 1))]
    if job.get("label"):
        parts += ["-l", shlex.quote(str(job["label"]))]
    if qmmm_parent is not None:
        parts += [qmmm_parent, "qmmm"]
        # Gaussian's QMMM settings are independent from the parent job
        # settings, so total state must appear in the nested QMMM scope. The
        # same explicit form is also valid and unambiguous for ORCA.
        if "charge_total" not in settings:
            subcommand_flags += [
                "--charge-total",
                str(job.get("charge", 0)),
            ]
        if "mult_total" not in settings:
            subcommand_flags += ["--mult-total", str(job.get("mult", 1))]
    else:
        parts.append(_subcommand(kind))
    parts += subcommand_flags
    return " ".join(parts)


def spec_to_commands(
    spec: dict[str, Any],
    default_project: str | None = None,
) -> list[str]:
    """Render a postprocessed workflow SPEC into ordered chemsmart commands."""
    if not isinstance(spec, dict) or spec.get("intent") != "workflow":
        return []
    commands: list[str] = []
    geom_of: dict[Any, str] = {}
    for job in spec.get("jobs", []):
        if not isinstance(job, dict):
            continue
        commands.append(_job_command(job, geom_of, default_project))
        geom_of[job.get("id")] = (
            str(job.get("label"))
            if job.get("label")
            else f"<{job.get('kind')}-output>"
        )
    return commands


def adapt(
    model_text: str | dict[str, Any],
    validate: bool = True,
    default_project: str | None = None,
) -> dict[str, Any]:
    """Parse, postprocess, render, and optionally validate a compact SPEC."""
    out: dict[str, Any] = {
        "intent": None,
        "spec": None,
        "commands": [],
        "valid": None,
        "errors": [],
        "message": None,
    }
    try:
        spec = (
            json.loads(model_text)
            if isinstance(model_text, str)
            else model_text
        )
    except Exception as exc:
        out["errors"].append(f"invalid JSON: {exc}")
        return out
    if not isinstance(spec, dict):
        out["errors"].append("SPEC must be a JSON object")
        return out

    spec = postprocess(spec)
    out["spec"] = spec
    out["intent"] = spec.get("intent")
    if spec.get("intent") != "workflow":
        out["message"] = spec.get("message")
        out["valid"] = True if validate else None
        return out

    try:
        out["commands"] = spec_to_commands(
            spec, default_project=default_project
        )
    except Exception as exc:
        out["valid"] = False if validate else None
        out["errors"].append(f"adapter render failed: {exc}")
        return out
    if validate:
        out["valid"], out["errors"] = _validate_all(out["commands"])
    return out


def _validate_all(commands: list[str]) -> tuple[bool, list[str]]:
    try:
        from chemsmart.cli.main import entry_point
    except Exception as exc:
        return False, [f"chemsmart CLI not importable: {exc}"]

    errors: list[str] = []
    for command in commands:
        try:
            tokens = shlex.split(command)
        except ValueError as exc:
            errors.append(f"{command!r}: {exc}")
            continue
        if not tokens or tokens[0] != "chemsmart":
            errors.append(f"{command!r}: must start with chemsmart")
            continue
        try:
            entry_point.make_context(
                "chemsmart",
                tokens[1:],
                resilient_parsing=True,
            )
        except Exception as exc:
            errors.append(f"{command!r}: {exc}")
            continue
        errors.extend(_semantic_errors(command, tokens))
    return (not errors), errors


def _semantic_errors(command: str, tokens: list[str]) -> list[str]:
    errors: list[str] = []
    try:
        from chemsmart.utils.utils import get_list_from_string_range
    except Exception:
        return errors
    for flag in _RANGE_FLAGS:
        if flag not in tokens:
            continue
        index = tokens.index(flag) + 1
        if index >= len(tokens):
            errors.append(f"{command!r}: {flag} requires a value")
            continue
        value = tokens[index]
        try:
            get_list_from_string_range(value)
        except Exception as exc:
            errors.append(
                f"{command!r}: {flag} failed runtime range parse: {exc}"
            )
    return errors


__all__ = ["adapt", "postprocess", "spec_to_commands"]
