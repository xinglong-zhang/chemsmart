"""chemsmart v8 adapter — bridges the local v8 spec-emission model to runnable chemsmart CLI commands.

Pipeline:  model text  ->  parse JSON spec  ->  postprocess (module-index cleanup)  ->  render chemsmart
command(s)  ->  validate against the real chemsmart CLI parser.

Drop-in for `chemsmart/agent/` (ships with `v8_kind_index.py` + `postprocess_v8.py`). The existing
SynthesisSession consumes a `{"command": "chemsmart ..."}` string; for the v8 model, call `adapt(text)`
and feed the resulting command(s) into that same validate/confirm/execute path.
"""
from __future__ import annotations
import json
import re

try:
    from . import v8_kind_index as KI
    from .postprocess_v8 import postprocess
except Exception:  # standalone use outside the package
    import v8_kind_index as KI
    from postprocess_v8 import postprocess

# kind suffix -> CLI subcommand. gaussian.sp/orca.sp -> singlepoint; tddft -> td; everything else == suffix.
# *.freq has NO standalone CLI subcommand -> route it as an `opt` job carrying a freq jobtype flag.
_SUBCMD = {"sp": "singlepoint", "tddft": "td"}
# v8 settings key -> CLI flag (names taken from chemsmart/cli/{gaussian,orca}/*.py @click.option)
_FLAG = {
    "additional_opt_options_in_route": "--additional-opt-options",
    "additional_route_parameters": "--additional-route-parameters",
    "direction": "--direction", "flat_irc": "--flat-irc", "inithess": "--inithess",
    "maxiter": "--maxiter", "hessmode": "--hessmode",
    "nstates": "--nstates", "states": "--states", "root": "--root", "eqsolv": "--eqsolv",
    "fragment_indices": "--fragment-indices",
    "num_confs_to_run": "--num-confs-to-run", "grouping_strategy": "--grouping-strategy",
    "num_groups": "--num-groups", "num_structures_to_run": "--num-structures-to-run",
    "proportion_structures_to_use": "--proportion-structures-to-use",
    "high_level_atoms": "--high-level-atoms", "low_level_atoms": "--low-level-atoms",
    "nimages": "--nimages", "joboption": "--joboption", "freeze_atoms": "--freeze-atoms",
    "recalc_hess": "--recalc-hess", "trust_radius": "--trust-radius", "tssearch_type": "--tssearch-type",
    "scan_definition": "--coordinates", "modred": "--coordinates",
}
import shlex


def _subcommand(kind: str) -> str:
    suf = kind.split(".", 1)[1]
    return _SUBCMD.get(suf, suf)


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, list):
        return shlex.quote(" ".join(str(x) for x in (sum(v, []) if v and isinstance(v[0], list) else v)))
    return shlex.quote(str(v))


def _job_command(job, geom_of):
    """Render one job to a chemsmart command (geom_of: id -> upstream output reference for chains)."""
    kind = job["kind"]; prog = kind.split(".", 1)[0]
    verb = "sub" if job.get("execution") == "submit" else "run"
    parts = ["chemsmart", verb]
    if job.get("server"):
        parts += ["-s", job["server"]]
    parts.append(prog)
    settings = dict(job.get("settings", {}) or {})
    freq_true = settings.pop("freq", None) is True
    # freq is a route token (verified B1: renders as " freq"); chemsmart has NO --freq flag and "opt_freq"
    # is not a real jobtype. Express opt+freq via --additional-route-parameters (combined with any existing).
    if freq_true and not kind.endswith(".freq"):
        aro = settings.get("additional_route_parameters")
        settings["additional_route_parameters"] = f"{aro} freq" if aro else "freq"
    for k, v in settings.items():
        flag = _FLAG.get(k)
        if flag:
            parts += [flag, _fmt(v)]
    # subcommand: *.freq has no standalone CLI subcommand -> run an opt job tagged as a frequency calc
    if kind.endswith(".freq"):
        parts += ["opt", "--jobtype", "freq"]
    else:
        parts.append(_subcommand(kind))
    src = job.get("file") or geom_of.get(job.get("geom_from"), "<upstream-geometry>")
    parts += ["-f", shlex.quote(src)]
    if kind == "orca.neb" and job.get("product_file"):
        parts += ["-e", shlex.quote(job["product_file"])]
    parts += ["-c", str(job["charge"]), "-m", str(job["mult"])]
    if job.get("label"):
        parts += ["-l", shlex.quote(job["label"])]
    return " ".join(parts)


def spec_to_commands(spec):
    """Postprocessed spec -> ordered list of chemsmart commands (one per job)."""
    if not isinstance(spec, dict) or spec.get("intent") != "workflow":
        return []
    geom_of = {}
    cmds = []
    for job in spec.get("jobs", []):
        cmds.append(_job_command(job, geom_of))
        # a downstream job's geometry comes from this job's optimized output (label-derived)
        geom_of[job.get("id")] = (job.get("label") or "<%s-output>" % job.get("kind"))
    return cmds


def adapt(model_text, validate=True):
    """Full bridge: parse -> postprocess -> render -> (optionally) validate against the chemsmart CLI.
    Returns {intent, spec, commands, valid, errors, message}."""
    out = {"intent": None, "spec": None, "commands": [], "valid": None, "errors": [], "message": None}
    try:
        spec = json.loads(model_text) if isinstance(model_text, str) else model_text
    except Exception as exc:
        out["errors"].append(f"invalid JSON: {exc}")
        return out
    spec = postprocess(spec)
    out["spec"] = spec
    out["intent"] = spec.get("intent")
    if spec.get("intent") != "workflow":
        out["message"] = spec.get("message")
        return out
    out["commands"] = spec_to_commands(spec)
    if validate:
        out["valid"], out["errors"] = _validate_all(out["commands"])
    return out


def _validate_all(commands):
    """Dry-parse each command against the real chemsmart CLI (chemsmart.cli.main.entry_point)."""
    try:
        from chemsmart.cli.main import entry_point
    except Exception as exc:  # chemsmart not importable in this env
        return None, [f"chemsmart CLI not importable: {exc}"]
    errs = []
    for cmd in commands:
        toks = shlex.split(cmd)
        if not toks or toks[0] != "chemsmart":
            errs.append(f"{cmd!r}: must start with chemsmart"); continue
        try:
            entry_point.make_context("chemsmart", toks[1:], resilient_parsing=True)
        except Exception as exc:
            errs.append(f"{cmd!r}: {exc}")
    return (not errs), errs


__all__ = ["postprocess", "spec_to_commands", "adapt"]
