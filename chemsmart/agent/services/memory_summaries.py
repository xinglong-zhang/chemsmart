"""Compact public-evidence summaries for conversation memory."""

from __future__ import annotations

import re
from collections import Counter
from typing import Any

_GAUSSIAN_ROUTE_RE = re.compile(r"^\s*#.*$", re.MULTILINE)
_ORCA_ROUTE_RE = re.compile(r"^\s*!.*$", re.MULTILINE)
_MOLECULE_REPR_RE = re.compile(r"Molecule<([^,>]+)")


def summarize_tool_result(
    payload: dict[str, Any],
    tool_call: dict[str, Any] | None,
) -> str | None:
    """Summarize a legacy tool-result payload."""

    tool = _string_value(payload.get("tool"))
    result = payload.get("payload")
    if tool is None or not isinstance(result, dict):
        return None
    args = tool_call.get("args") if isinstance(tool_call, dict) else {}
    if not isinstance(args, dict):
        args = {}
    if tool == "build_molecule":
        return _build_molecule_summary(args, result)
    if tool == "recommend_method":
        method = _method_summary(result)
        return f"recommend_method suggested {method}." if method else None
    if tool in _PROJECT_YAML_TOOLS:
        return _project_yaml_summary(tool, result)
    if tool in {"build_gaussian_settings", "build_orca_settings"}:
        method = _method_summary(result)
        return f"{tool} prepared {method} settings." if method else None
    if tool == "build_job":
        return _build_job_summary(args, result)
    if tool == "dry_run_input":
        return _dry_run_summary(result)
    if tool in {"execute_chemsmart_command", "inspect_calculation"}:
        return _calculation_summary(tool, result)
    if tool == "extract_optimized_geometry":
        return _geometry_summary(_molecule_formula(result))
    if tool == "validate_runtime":
        status = _string_value(result.get("ok"))
        return (
            f"validate_runtime reported {status}."
            if status and status != "ok"
            else None
        )
    return None


def summarize_tool_use_result(
    payload: dict[str, Any],
    req: dict[str, Any] | None,
) -> str | None:
    """Summarize a run-loop tool-use result in display-result form."""

    if not isinstance(payload, dict):
        return None
    tool = _string_value(payload.get("tool"))
    if tool is None or payload.get("status") not in ("ok", "partial"):
        return None
    inner = payload.get("payload")
    args = req.get("args") if isinstance(req, dict) else {}
    if not isinstance(args, dict):
        args = {}
    if tool == "build_molecule":
        return _display_molecule_summary(args, inner)
    if tool == "recommend_method":
        method = _method_summary(inner) if isinstance(inner, dict) else None
        return f"recommend_method suggested {method}." if method else None
    if tool in _PROJECT_YAML_TOOLS:
        return _project_yaml_summary(
            tool,
            inner if isinstance(inner, dict) else {},
        )
    if tool in {"build_gaussian_settings", "build_orca_settings"}:
        method = _method_summary(
            {
                "functional": args.get("functional"),
                "basis": args.get("basis"),
                "ab_initio": args.get("ab_initio"),
                "solvent_id": args.get("solvent_id"),
            }
        )
        return f"{tool} prepared {method} settings." if method else None
    if tool == "build_job":
        return _display_job_summary(args, inner)
    if tool == "dry_run_input":
        return _display_dry_run_summary(inner)
    if tool == "extract_optimized_geometry":
        return _geometry_summary(_display_molecule_formula(inner))
    if tool == "validate_runtime":
        return _display_runtime_summary(inner)
    if tool == "ssh_probe":
        return _ssh_probe_summary(args, inner, payload.get("status"))
    if tool == "scheduler_query":
        return _scheduler_summary(args, inner)
    if tool == "log_tail":
        return _log_tail_summary(args, inner)
    if tool == "read":
        return _read_summary(args, inner)
    return None


def summarize_ask_user(payload: dict[str, Any]) -> str | None:
    question = _string_value(payload.get("question"))
    if question is None:
        return None
    options = payload.get("options")
    option_list = (
        [
            option.strip()
            for option in options
            if isinstance(option, str) and option.strip()
        ]
        if isinstance(options, list)
        else []
    )
    if option_list:
        return (
            f"ask_user requested clarification: {question} "
            f"Options: {', '.join(option_list[:4])}."
        )
    return f"ask_user requested clarification: {question}."


def summarize_ask_user_answer(payload: dict[str, Any]) -> str | None:
    answer = _string_value(payload.get("answer"))
    return f"User answered clarification: {answer}." if answer else None


_PROJECT_YAML_TOOLS = {
    "extract_project_protocol",
    "render_project_yaml",
    "validate_project_yaml",
    "critic_project_yaml",
    "write_project_yaml",
}


def _build_molecule_summary(
    args: dict[str, Any], payload: dict[str, Any]
) -> str:
    source = _string_value(args.get("filepath")) or _string_value(
        args.get("smiles")
    )
    formula = _molecule_formula(payload)
    if source and formula:
        return f"build_molecule loaded {formula} from {source}."
    if source:
        return f"build_molecule loaded the source from {source}."
    if formula:
        return f"build_molecule produced molecule {formula}."
    return "build_molecule produced a reusable molecule."


def _build_job_summary(args: dict[str, Any], payload: dict[str, Any]) -> str:
    kind = _string_value(args.get("kind"))
    label = _string_value(payload.get("label")) or _string_value(
        args.get("label")
    )
    pieces = ["build_job prepared"]
    if kind:
        pieces.append(kind)
    if label:
        pieces.append(f"label={label}")
    return " ".join(pieces) + "."


def _dry_run_summary(payload: dict[str, Any]) -> str:
    command = _string_value(payload.get("command"))
    inputfile = _string_value(payload.get("inputfile"))
    route = _extract_route_line(payload.get("content"))
    pieces = ["dry_run_input wrote"]
    if inputfile:
        pieces.append(inputfile)
    if command:
        pieces.append(f"for command {command}")
    if route:
        pieces.append(f"with route {route}")
    return " ".join(pieces) + "."


def _calculation_summary(tool: str, payload: dict[str, Any]) -> str | None:
    calculation = payload.get("calculation")
    if not isinstance(calculation, dict):
        return None
    run_id = _string_value(calculation.get("run_id")) or "latest"
    status = _string_value(calculation.get("status")) or "unknown"
    output_path = _string_value(calculation.get("output_path"))
    energy = calculation.get("energy")
    pieces = [f"{tool} recorded {run_id} status={status}"]
    if output_path:
        pieces.append(f"output={output_path}")
    if isinstance(energy, (int, float)):
        pieces.append(f"energy={float(energy):.12f} Eh")
    return ", ".join(pieces) + "."


def _display_molecule_summary(args: dict[str, Any], inner: Any) -> str:
    source = _string_value(args.get("filepath")) or _string_value(
        args.get("smiles")
    )
    formula = _display_molecule_formula(inner)
    if source and formula:
        return f"build_molecule loaded {formula} from {source}."
    if source:
        return f"build_molecule loaded the source from {source}."
    if formula:
        return f"build_molecule produced molecule {formula}."
    return "build_molecule produced a reusable molecule."


def _display_molecule_formula(inner: Any) -> str | None:
    if not isinstance(inner, dict):
        return None
    summary = inner.get("summary")
    text = _string_value(
        summary.get("repr") if isinstance(summary, dict) else None
    )
    if not text:
        return None
    match = _MOLECULE_REPR_RE.search(text)
    return match.group(1) if match else None


def _display_job_summary(args: dict[str, Any], inner: Any) -> str:
    kind = _string_value(args.get("kind"))
    label = _string_value(args.get("label"))
    if not label and isinstance(inner, dict):
        summary = inner.get("summary")
        text = _string_value(
            summary.get("repr") if isinstance(summary, dict) else None
        )
        if text and (match := re.search(r"label=([^,>]+)", text)):
            label = match.group(1).strip()
    pieces = ["build_job prepared"]
    if kind:
        pieces.append(kind)
    if label:
        pieces.append(f"label={label}")
    return " ".join(pieces) + "."


def _display_dry_run_summary(inner: Any) -> str:
    summary = inner.get("summary") if isinstance(inner, dict) else {}
    summary = summary if isinstance(summary, dict) else {}
    return _dry_run_summary(summary)


def _display_runtime_summary(inner: Any) -> str | None:
    ok_value = None
    if isinstance(inner, dict):
        summary = inner.get("summary")
        ok_value = _string_value(
            summary.get("ok") if isinstance(summary, dict) else inner.get("ok")
        )
    if ok_value and ok_value != "ok":
        return f"validate_runtime reported {ok_value}."
    return None


def _ssh_probe_summary(args: dict[str, Any], inner: Any, status: Any) -> str:
    payload = inner if isinstance(inner, dict) else {}
    server = _string_value(args.get("server"))
    probe = _string_value(args.get("probe_name")) or _string_value(
        payload.get("probe")
    )
    observed = _string_value(payload.get("scheduler")) or _string_value(status)
    pieces = ["ssh_probe"]
    if probe:
        pieces.append(probe)
    if server:
        pieces.append(f"on {server}")
    text = " ".join(pieces)
    return f"{text} → {observed}" if observed else text


def _scheduler_summary(args: dict[str, Any], inner: Any) -> str:
    payload = inner if isinstance(inner, dict) else {}
    server = _string_value(args.get("server"))
    pieces = ["scheduler_query"]
    if server:
        pieces.append(f"on {server}:")
    else:
        pieces[-1] += ":"
    details = []
    job_id = _string_value(payload.get("job_id")) or _string_value(
        args.get("job_id")
    )
    state = _string_value(payload.get("state"))
    queue = _string_value(payload.get("queue")) or _string_value(
        payload.get("partition_or_queue")
    )
    if job_id:
        details.append(f"job {job_id}")
    if state:
        details.append(f"state={state}")
    if queue:
        details.append(f"queue={queue}")
    return " ".join([*pieces, *details]).rstrip(":")


def _log_tail_summary(args: dict[str, Any], inner: Any) -> str:
    payload = inner if isinstance(inner, dict) else {}
    path = _string_value(args.get("path")) or _string_value(
        payload.get("path")
    )
    lines_returned = _coerce_int(payload.get("lines_returned"))
    errors = payload.get("errors")
    error_count = len(errors) if isinstance(errors, list) else 0
    top_kind = None
    if isinstance(errors, list) and errors and isinstance(errors[0], dict):
        top_kind = _string_value(errors[0].get("kind"))
    details = []
    if lines_returned is not None:
        details.append(f"{lines_returned}L")
    details.append(f"{error_count} errors")
    if top_kind:
        details[-1] += f": {top_kind}"
    text = "log_tail" + (f" {path}" if path else "")
    return text + (f" ({', '.join(details)})" if details else "")


def _read_summary(args: dict[str, Any], inner: Any) -> str:
    payload = inner if isinstance(inner, dict) else {}
    path = _string_value(args.get("path")) or _string_value(
        payload.get("path")
    )
    start = _coerce_int(payload.get("start_line"))
    end = _coerce_int(payload.get("end_line"))
    total = _coerce_int(payload.get("total_lines"))
    text = "read" + (f" {path}" if path else "")
    if start is not None and end is not None and total is not None:
        text += f" L{start}-{end}/{total}"
    return text


def _molecule_formula(payload: dict[str, Any]) -> str | None:
    symbols = payload.get("symbols")
    if not isinstance(symbols, list) or not symbols:
        return None
    counts = Counter(
        symbol
        for symbol in symbols
        if isinstance(symbol, str) and symbol.strip()
    )
    if not counts:
        return None
    return "".join(
        symbol if count == 1 else f"{symbol}{count}"
        for symbol, count in sorted(counts.items())
    )


def _method_summary(payload: dict[str, Any]) -> str | None:
    functional = _string_value(payload.get("functional"))
    basis = _string_value(payload.get("basis"))
    ab_initio = _string_value(payload.get("ab_initio"))
    solvent = _string_value(payload.get("solvent_id"))
    method = None
    if ab_initio and basis:
        method = f"{ab_initio}/{basis}"
    elif functional and basis:
        method = f"{functional}/{basis}"
    elif functional:
        method = functional
    return f"{method} in {solvent}" if method and solvent else method


def _project_yaml_summary(tool: str, payload: dict[str, Any]) -> str | None:
    project = _string_value(payload.get("project_name")) or "project"
    program = _string_value(payload.get("program")) or "program"
    verdict = _string_value(payload.get("verdict"))
    if verdict is None and isinstance(payload.get("validation"), dict):
        verdict = _string_value(payload["validation"].get("verdict"))
    if tool == "extract_project_protocol":
        return f"extract_project_protocol extracted {program}:{project} method facts."
    if tool == "render_project_yaml":
        return f"render_project_yaml produced a {program}:{project} YAML candidate."
    if tool == "write_project_yaml":
        path = _string_value(payload.get("written_path"))
        if path:
            return f"write_project_yaml wrote {program}:{project} to {path}."
        return f"write_project_yaml did not write {program}:{project}."
    return (
        f"{tool} returned {verdict} for {program}:{project}."
        if verdict
        else None
    )


def _geometry_summary(formula: str | None) -> str:
    if formula:
        return (
            "extract_optimized_geometry recovered optimized geometry "
            f"for {formula}."
        )
    return "extract_optimized_geometry recovered optimized geometry."


def _extract_route_line(content: Any) -> str | None:
    if not isinstance(content, str):
        return None
    match = _GAUSSIAN_ROUTE_RE.search(content) or _ORCA_ROUTE_RE.search(
        content
    )
    return _truncate(match.group(0).strip(), 120) if match else None


def _string_value(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _truncate(text: str, limit: int) -> str:
    if len(text) <= limit:
        return text
    return text[: max(0, limit - 1)].rstrip() + "…"


__all__ = [
    "summarize_ask_user",
    "summarize_ask_user_answer",
    "summarize_tool_result",
    "summarize_tool_use_result",
]
