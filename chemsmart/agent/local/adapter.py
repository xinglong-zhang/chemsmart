"""Translate local planner JSON to ``chemsmart agent ask`` synthesis schema.

The v13.1 model emits compact ``{"intent":"workflow","jobs":[...]}`` SPECs.
Those route through :mod:`chemsmart.agent.v8_adapter`, which postprocesses and
validates against the real CLI/parser contract. The older V4 ``steps`` shape is
kept as a compatibility fallback for local experiments and old tests.

The conversion is conservative: when arguments are missing or the kind is not
yet wired into the CLI, the adapter returns ``status="needs_clarification"``
with the missing slots listed in ``missing_info`` so the caller can re-prompt.
"""

from __future__ import annotations

import shlex
from typing import Any

# build_job.kind -> (chemsmart program subcommand, job subcommand)
_KIND_TO_CLI: dict[str, tuple[str, str]] = {
    "gaussian.sp": ("gaussian", "sp"),
    "gaussian.opt": ("gaussian", "opt"),
    "gaussian.ts": ("gaussian", "ts"),
    "gaussian.irc": ("gaussian", "irc"),
    "gaussian.scan": ("gaussian", "modred"),
    "gaussian.modred": ("gaussian", "modred"),
    "gaussian.tddft": ("gaussian", "td"),
    "orca.sp": ("orca", "sp"),
    "orca.opt": ("orca", "opt"),
    "orca.ts": ("orca", "ts"),
    "orca.irc": ("orca", "irc"),
    "orca.scan": ("orca", "scan"),
}

# Kinds the V4 LoRA can emit but the CLI does not (yet) accept as drop-in
# top-level subcommands. These force a needs_clarification with a helpful
# message rather than an invalid command.
_UNSUPPORTED_KINDS: frozenset[str] = frozenset(
    {
        "gaussian.resp",
        "gaussian.nci",
        "gaussian.dias",
        "gaussian.crest",
        "gaussian.traj",
        "gaussian.qmmm",
        "gaussian.qrc",
        "gaussian.wbi",
        "gaussian.freq",
        "orca.modred",
        "orca.neb",
        "orca.qmmm",
        "orca.qrc",
    }
)


def plan_to_synthesis_result(
    plan: dict[str, Any],
    user_query: str,
    use_submit: bool | None = None,
    default_project: str | None = None,
) -> dict[str, Any]:
    """Convert local planner JSON into a ``SynthesisSession``-shaped result.

    Args:
        plan: Parsed planner dict, ideally already post-processed.
        user_query: Original user query (used only for the explanation).
        use_submit: Force ``sub`` (True) or ``run`` (False). When ``None`` the
            decision is taken from the presence of a ``submit_hpc`` step.
        default_project: Runtime-owned chemsmart project name to inject into
            compact-SPEC commands. The model must not emit this field.

    Returns:
        A dict matching the contract of
        :func:`chemsmart.agent.synthesis._normalize_result`. ``status`` is
        ``"ready"`` only when a single legal chemsmart command can be emitted.
    """
    if _is_compact_spec(plan):
        return _compact_spec_to_synthesis_result(
            plan,
            default_project=default_project,
        )

    intent = str(plan.get("intent") or "").lower()
    if intent in {"decline", "infeasible", "out_of_scope"}:
        return _infeasible(
            plan.get("rationale") or "Planner declined the request."
        )
    if intent in {"advisory", "needs_clarification", "clarify"}:
        return _needs_clarification(
            plan.get("rationale") or "Planner needs more information.",
            _extract_missing_info(plan),
        )

    steps = plan.get("steps")
    if not isinstance(steps, list) or not steps:
        return _needs_clarification(
            "Planner produced no executable steps.",
            ["intent", "steps"],
        )

    job_step = _first_build_job(steps)
    if job_step is None:
        return _needs_clarification(
            "Planner did not emit a build_job step; cannot synthesize a CLI "
            "command.",
            ["build_job.kind"],
        )

    job_args = job_step.get("args") or {}
    kind = str(job_args.get("kind", "")).strip().lower()
    if kind in _UNSUPPORTED_KINDS:
        return _needs_clarification(
            f"Planner emitted kind={kind!r}, which has no direct chemsmart "
            "CLI subcommand. Re-phrase the request as one of: "
            f"{sorted(_KIND_TO_CLI)}.",
            [f"kind:{kind}->cli"],
        )
    if kind not in _KIND_TO_CLI:
        return _needs_clarification(
            f"Unknown build_job.kind {kind!r}.",
            ["build_job.kind"],
        )

    program, subcommand = _KIND_TO_CLI[kind]
    filepath = _extract_filepath(steps)
    if not filepath:
        return _needs_clarification(
            "Planner did not specify an input molecule filepath.",
            ["molecule.filepath"],
        )

    settings = _extract_settings(steps)
    submit_step = _find_submit_step(steps)
    submitting = (
        use_submit if use_submit is not None else submit_step is not None
    )
    tokens: list[str] = ["chemsmart", "sub" if submitting else "run"]
    if submitting and submit_step is not None:
        server = str((submit_step.get("args") or {}).get("server", "")).strip()
        if server:
            tokens += ["-s", server]
        cores = (submit_step.get("args") or {}).get("ncpus")
        if isinstance(cores, int) and cores > 0:
            tokens += ["-n", str(cores)]
    tokens.append(program)
    tokens += _settings_flags(program, settings)
    tokens += ["-f", filepath]
    label = str(job_args.get("label", "")).strip()
    if label:
        tokens += ["-l", label]
    tokens.append(subcommand)

    command = shlex.join(tokens)
    explanation = _build_explanation(plan, user_query, kind, settings)
    return {
        "status": "ready",
        "command": command,
        "explanation": explanation,
        "confidence": _confidence_for(plan, settings),
        "missing_info": [],
        "alternatives": [],
    }


def _is_compact_spec(plan: dict[str, Any]) -> bool:
    intent = plan.get("intent")
    return isinstance(intent, str) and (
        intent in {"workflow", "advisory", "decline", "chitchat"}
    ) and ("jobs" in plan or "message" in plan)


def _compact_spec_to_synthesis_result(
    plan: dict[str, Any],
    default_project: str | None = None,
) -> dict[str, Any]:
    from chemsmart.agent.v8_adapter import adapt

    adapted = adapt(plan, validate=True, default_project=default_project)
    intent = str(adapted.get("intent") or plan.get("intent") or "")
    if intent != "workflow":
        message = str(adapted.get("message") or plan.get("message") or "")
        status = "infeasible" if intent == "decline" else "needs_clarification"
        if intent == "chitchat":
            status = "infeasible"
        return {
            "status": status,
            "command": "",
            "explanation": message or "No executable workflow was requested.",
            "confidence": "high",
            "missing_info": [],
            "alternatives": [],
        }

    commands = [
        command
        for command in adapted.get("commands", [])
        if isinstance(command, str) and command.strip()
    ]
    if not commands:
        errors = adapted.get("errors") or ["no commands rendered"]
        return _needs_clarification(
            "Planner SPEC could not be rendered: " + "; ".join(map(str, errors)),
            ["jobs"],
        )
    if adapted.get("valid") is False:
        return _needs_clarification(
            "Rendered command failed chemsmart validation: "
            + "; ".join(map(str, adapted.get("errors") or [])),
            ["validated_cli_command"],
        )

    explanation = "Prepared chemsmart command from compact SPEC."
    if len(commands) > 1:
        explanation = (
            "Prepared a multi-step chemsmart workflow from compact SPEC."
        )
    return {
        "status": "ready",
        "command": commands[0],
        "explanation": explanation,
        "confidence": "high" if adapted.get("valid") else "medium",
        "missing_info": [],
        "alternatives": commands[1:],
    }


def _first_build_job(steps: list[Any]) -> dict[str, Any] | None:
    for step in steps:
        if isinstance(step, dict) and step.get("tool") == "build_job":
            return step
    return None


def _find_submit_step(steps: list[Any]) -> dict[str, Any] | None:
    for step in steps:
        if isinstance(step, dict) and step.get("tool") == "submit_hpc":
            return step
    return None


def _extract_filepath(steps: list[Any]) -> str:
    for step in steps:
        if not isinstance(step, dict) or step.get("tool") != "build_molecule":
            continue
        args = step.get("args") or {}
        filepath = args.get("filepath")
        if isinstance(filepath, str) and filepath.strip():
            return filepath.strip()
    return ""


def _extract_settings(steps: list[Any]) -> dict[str, Any]:
    for step in steps:
        if not isinstance(step, dict):
            continue
        tool = step.get("tool", "")
        if tool in {"build_gaussian_settings", "build_orca_settings"}:
            args = step.get("args") or {}
            if isinstance(args, dict):
                return args
    return {}


def _settings_flags(program: str, settings: dict[str, Any]) -> list[str]:
    flags: list[str] = []
    functional = settings.get("functional")
    if isinstance(functional, str) and functional:
        flags += ["-x", functional]
    basis = settings.get("basis")
    if isinstance(basis, str) and basis:
        flags += ["-b", basis]
    charge = settings.get("charge")
    if isinstance(charge, int):
        flags += ["-c", str(charge)]
    multiplicity = settings.get("multiplicity")
    if isinstance(multiplicity, int):
        flags += ["-m", str(multiplicity)]
    solvent_model = settings.get("solvent_model")
    solvent_id = settings.get("solvent_id")
    if isinstance(solvent_model, str) and solvent_model:
        flags += ["--solvent-model", solvent_model]
    if isinstance(solvent_id, str) and solvent_id:
        flags += ["--solvent-id", solvent_id]
    return flags


def _extract_missing_info(plan: dict[str, Any]) -> list[str]:
    missing = plan.get("missing_info")
    if isinstance(missing, list):
        return [str(item) for item in missing if isinstance(item, str)]
    return []


def _build_explanation(
    plan: dict[str, Any],
    user_query: str,
    kind: str,
    settings: dict[str, Any],
) -> str:
    pieces: list[str] = []
    rationale = plan.get("rationale")
    if isinstance(rationale, str) and rationale.strip():
        pieces.append(rationale.strip())
    method = settings.get("functional")
    basis = settings.get("basis")
    if method and basis:
        pieces.append(f"Method: {method}/{basis} ({kind}).")
    if user_query.strip():
        pieces.append(f"From: {user_query.strip()}")
    return " ".join(pieces)[:600]


def _confidence_for(plan: dict[str, Any], settings: dict[str, Any]) -> str:
    if not settings.get("functional") or not settings.get("basis"):
        return "low"
    intent = str(plan.get("intent") or "").lower()
    if intent in {"workflow", "execute"}:
        return "high"
    return "medium"


def _needs_clarification(reason: str, missing: list[str]) -> dict[str, Any]:
    return {
        "status": "needs_clarification",
        "command": "",
        "explanation": reason,
        "confidence": "low",
        "missing_info": missing,
        "alternatives": [],
    }


def _infeasible(reason: str) -> dict[str, Any]:
    return {
        "status": "infeasible",
        "command": "",
        "explanation": reason,
        "confidence": "low",
        "missing_info": [],
        "alternatives": [],
    }
