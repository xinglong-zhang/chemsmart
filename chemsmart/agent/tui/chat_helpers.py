"""Pure projections and compatibility helpers for the chat screen."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from chemsmart.agent.provider_config import (
    AgentProviderConfig,
    AgentProviderConfigError,
    load_active_provider_config,
)
from chemsmart.agent.tui.chat_models import ReadyCommand
from chemsmart.settings.workspace_project import WorkspaceProjectStatus


class _ExecuteOverrideRegistry:
    """Wraps ToolRegistry to inject execute=True + transport into submit_hpc.

    Used by execute_agent_session so the LLM does not need to know about the
    internal execute flag — the framework injects it transparently.
    """

    def __init__(self, base, transport) -> None:
        self._base = base
        self._transport = transport

    def __getattr__(self, name: str):
        return getattr(self._base, name)

    def call(self, tool_name: str, args: dict) -> object:
        if tool_name == "submit_hpc":
            args = {**args, "execute": True, "transport": self._transport}
        return self._base.call(tool_name, args)


def _extract_execute_context(session_dir: Path) -> dict[str, str | None]:
    """Read the session decision log and pull out the job handle and inputfile.

    Returns a dict with 'job_handle_id', 'inputfile', and 'server_name'.
    All values may be None if not found.
    """
    log_path = session_dir / "decision_log.jsonl"
    job_handle_id: str | None = None
    inputfile: str | None = None
    server_name: str | None = None
    if not log_path.exists():
        return {
            "job_handle_id": job_handle_id,
            "inputfile": inputfile,
            "server_name": server_name,
        }
    try:
        for raw in log_path.read_text(encoding="utf-8").splitlines():
            if not raw.strip():
                continue
            entry = json.loads(raw)
            kind = entry.get("kind")
            payload = entry.get("payload") or {}
            if kind == "tool_use_result" and isinstance(payload, dict):
                tool = payload.get("tool")
                if tool == "build_job" and payload.get("status") == "ok":
                    job_handle_id = str(payload["handle_id"])
                elif tool == "dry_run_input" and payload.get("status") == "ok":
                    inner = payload.get("payload") or {}
                    summary = inner.get("summary") or {}
                    inputfile = summary.get("inputfile") or None
            elif kind == "tool_use_request" and isinstance(payload, dict):
                if payload.get("tool") == "submit_hpc":
                    args = payload.get("args") or {}
                    server_name = args.get("server_name") or server_name
    except Exception:
        pass
    return {
        "job_handle_id": job_handle_id,
        "inputfile": inputfile,
        "server_name": server_name,
    }


def _latest_project_yaml_candidate(
    session_dir: Path | None,
) -> dict[str, object] | None:
    if session_dir is None:
        return None
    log_path = session_dir / "decision_log.jsonl"
    if not log_path.exists():
        return None

    candidate: dict[str, object] | None = None
    candidate_yaml: object | None = None
    validated = False
    validation_requests: dict[str, dict[str, object]] = {}
    try:
        for raw in log_path.read_text(encoding="utf-8").splitlines():
            if not raw.strip():
                continue
            entry = json.loads(raw)
            payload = entry.get("payload") or {}
            if not isinstance(payload, dict):
                continue
            tool = payload.get("tool")
            call_id = str(payload.get("provider_call_id") or "")
            if (
                entry.get("kind") == "tool_use_request"
                and tool == "validate_project_yaml"
                and call_id
            ):
                args = payload.get("args") or {}
                yaml_text = args.get("yaml_text") if isinstance(args, dict) else None
                if isinstance(yaml_text, str) and not yaml_text.strip():
                    continue
                if not isinstance(yaml_text, (str, dict)):
                    continue
                validation_requests[call_id] = {
                    "project_name": str(
                        args.get("project_name") or "project"
                    ),
                    "program": str(args.get("program") or "gaussian"),
                    "yaml_text": yaml_text,
                }
                continue
            if payload.get("status") not in {None, "ok"}:
                continue

            result = payload.get("payload") or {}
            if isinstance(result, dict) and isinstance(
                result.get("summary"), dict
            ):
                result = result["summary"]
            if not isinstance(result, dict):
                continue
            if tool == "validate_project_yaml":
                if result.get("verdict") not in {"ok", "warn", "reject"}:
                    continue
                request_candidate = validation_requests.pop(call_id, None)
                validated = result.get("verdict") in {"ok", "warn"}
                if validated and request_candidate is not None:
                    candidate = request_candidate
                    candidate_yaml = request_candidate["yaml_text"]
                    continue
                # Older logs recorded only results. They may validate a prior
                # accepted render, but never a rejected or missing candidate.
                if candidate is None:
                    validated = False
                    continue
                if not validated:
                    candidate = None
                    candidate_yaml = None
                continue
            if tool != "render_project_yaml":
                continue
            yaml_text = result.get("yaml_text")
            if not isinstance(yaml_text, str) or not yaml_text.strip():
                continue
            validation = result.get("validation")
            render_rejected = result.get("ok") is False or (
                isinstance(validation, dict)
                and validation.get("verdict") == "reject"
            )
            if render_rejected:
                candidate = None
                candidate_yaml = None
                validated = False
                continue
            # A re-render of the identical candidate (build-mode over-iteration)
            # must NOT discard a prior successful validation — only a genuinely
            # new/changed candidate resets the validated flag and needs
            # re-validation before it can be written.
            if yaml_text == candidate_yaml:
                continue
            candidate = {
                "project_name": str(result.get("project_name") or "project"),
                "program": str(result.get("program") or "gaussian"),
                "yaml_text": yaml_text,
            }
            candidate_yaml = yaml_text
            validated = False
    except Exception:
        return None
    return candidate if validated else None


def _yaml_footer_label(status: WorkspaceProjectStatus) -> str:
    if status.loaded:
        return f"YAML OK {status.program}:{status.project}"
    if status.candidates:
        return f"YAML SELECT {len(status.candidates)}"
    return "YAML MISSING"


def _find_project_yaml_candidate_for_write(
    session_root: Path,
    *,
    preferred_session_dir: Path | None = None,
) -> dict[str, object] | None:
    candidate = _latest_project_yaml_candidate(preferred_session_dir)
    if candidate is not None:
        return candidate
    if not session_root.exists():
        return None
    try:
        session_dirs = sorted(
            [
                path
                for path in session_root.iterdir()
                if path.is_dir() and not path.name.startswith(".")
            ],
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
    except OSError:
        return None
    for session_dir in session_dirs:
        if (
            preferred_session_dir is not None
            and session_dir == preferred_session_dir
        ):
            continue
        candidate = _latest_project_yaml_candidate(session_dir)
        if candidate is not None:
            return candidate
    return None


def _tool_use_payload_summary(payload: object) -> str | None:
    if not isinstance(payload, dict):
        return None
    if payload.get("ok") is False:
        error = payload.get("error")
        if isinstance(error, dict):
            return str(error.get("message") or "") or None
    if "handle_id" in payload:
        return f"returned {payload['handle_id']}"
    if "inputfile" in payload:
        return str(payload.get("inputfile"))
    if "ok" in payload and payload.get("ok") is True:
        return "completed successfully"
    return None


def _tool_use_summary_payload(payload: object) -> dict[str, object] | None:
    if not isinstance(payload, dict):
        return None
    summary = payload.get("summary")
    if isinstance(summary, dict):
        return summary
    if "content" in payload or "inputfile" in payload:
        return payload
    return None


def _public_tool_result_payload(payload: object) -> dict[str, object] | None:
    """Project tool output without raw model responses or private reasoning."""

    if not isinstance(payload, dict):
        return None
    summary = payload.get("summary")
    if isinstance(summary, dict):
        payload = summary
    allowed = {
        "ok",
        "status",
        "error",
        "message",
        "path",
        "handle_id",
        "inputfile",
        "command",
        "semantic",
        "verdict",
        "failed_rule_ids",
        "issues",
        "job_id",
        "server",
        "returncode",
        "cli_grounded",
        "cli_grounding_issue",
        "calculation",
    }
    projected = {
        key: value for key, value in payload.items() if key in allowed
    }
    return projected or None


def _jobs_sort_key(job: dict) -> tuple[int, str, str]:
    order = {"running": 0, "queued": 1, "failed": 2, "cancelled": 3, "done": 4}
    return (
        order.get(str(job.get("status")), 9),
        str(job.get("raw_started") or ""),
        str(job.get("job_id") or ""),
    )


def _load_tui_provider_config() -> AgentProviderConfig | None:
    # Mirror ``get_provider``: the active ``agent.yaml`` config is the source of
    # truth for which provider the harness/synthesis path actually builds. A
    # stale ``AI_PROVIDER`` env var must NOT shadow it — otherwise the TUI picks
    # a mode that does not match the provider it will instantiate (e.g. defaults
    # to run/harness while ``get_provider`` returns a local provider, crashing
    # with "Unsupported provider 'local'").
    try:
        return load_active_provider_config()
    except AgentProviderConfigError:
        return None


def _provider_type_label(config: AgentProviderConfig | None) -> str:
    if config is None:
        return "offline"
    return str(config.type or "offline").strip().lower() or "offline"


def _provider_model_label(config: AgentProviderConfig | None) -> str:
    if config is None:
        return "auto"
    return str(config.model or "auto").strip() or "auto"


def _new_synthesis_artifact_dir(
    session_root: Path,
    *,
    provider_type: str,
) -> Path:
    lane = "local" if provider_type == "local" else "api"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = session_root / lane / "synthesis" / stamp
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_synthesis_artifact(
    artifact_dir: Path,
    payload: dict[str, object],
) -> None:
    try:
        artifact_dir.mkdir(parents=True, exist_ok=True)
        (artifact_dir / "synthesis_turn.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True, default=str),
            encoding="utf-8",
        )
    except OSError:
        return


def _default_interaction_mode(config: AgentProviderConfig | None) -> str:
    return (
        "synthesis"
        if config is not None and config.type == "local"
        else "unified"
    )


def _ready_command_hint(ready: ReadyCommand | None) -> str:
    if ready is None:
        return "Command shown, but execution evidence is incomplete"
    if ready.action == "run":
        return "Command validated — /run to execute locally"
    return "Command validated — /submit to submit for real"


def _is_calculation_diagnostic_request(request: str) -> bool:
    text = str(request or "").lower()
    return any(
        marker in text
        for marker in (
            "diagnose the result",
            "inspect the result",
            "analyze the result",
            "calculation result",
            "check the output",
            "계산 결과",
            "결과를 진단",
            "결과 분석",
            "출력 파일 확인",
            "诊断结果",
            "分析计算结果",
        )
    )


def _calculation_diagnostic_summary(calculation: dict[str, object]) -> str:
    program = str(calculation.get("program") or "calculation").upper()
    kind = str(calculation.get("kind") or "job").upper()
    status = str(calculation.get("status") or "parsed")
    lines = [f"{program} {kind} status: `{status}`"]
    energy = calculation.get("energy")
    if isinstance(energy, (int, float)):
        lines.append(f"- Final electronic energy: `{float(energy):.12f} Eh`")
    cycles = calculation.get("scf_cycles")
    if isinstance(cycles, int):
        lines.append(f"- SCF convergence: `{cycles} cycles`")
    imag = list(calculation.get("imag_freqs") or [])
    if imag:
        lines.append(
            "- Imaginary frequencies: `"
            + ", ".join(f"{float(value):.1f}" for value in imag)
            + " cm^-1`"
        )
    normal = calculation.get("normal_termination")
    if normal is not None:
        lines.append(
            "- Program termination: `normal`"
            if normal
            else "- Program termination: `not confirmed`"
        )
    output_path = str(calculation.get("output_path") or "")
    if output_path:
        lines.append(f"- Output: `{output_path}`")
    if status != "completed":
        error = str(calculation.get("error") or "").strip()
        if error:
            lines.extend(["", "Relevant diagnostic context:", "```text"])
            lines.extend(error.splitlines()[-40:])
            lines.append("```")
    return "\n".join(lines)


def _decision_trace_dict(
    synthesis: dict[str, object],
) -> dict[str, object] | None:
    trace = synthesis.get("decision_trace")
    return trace if isinstance(trace, dict) and trace else None


def _final_command_text(
    *,
    command: str,
) -> str:
    return "\n".join(("```bash", command, "```"))


def _command_details_text(
    *,
    explanation: str,
    confidence: str,
    project: str,
) -> str:
    parts = [f"confidence: `{confidence}`"]
    if project:
        parts.extend(["", f"active project: `{project}`"])
    if explanation:
        parts.extend(["", explanation])
    return "\n".join(parts)


def _final_answer_text(*, command: str, explanation: str) -> str:
    parts: list[str] = []
    if command:
        parts.extend(["Analyzed command:", "", "```bash", command, "```", ""])
    parts.append(explanation)
    return "\n".join(parts)


def _format_semantic_result(semantic: dict[str, object] | None) -> str:
    if not semantic:
        return ""
    verdict = str(semantic.get("verdict") or "unknown")
    failed = semantic.get("failed_rule_ids") or []
    if isinstance(failed, list) and failed:
        failed_text = ", ".join(str(item) for item in failed)
    else:
        failed_text = "none"
    issues = semantic.get("issues") or []
    issue_lines = []
    if isinstance(issues, list):
        for issue in issues:
            if not isinstance(issue, dict):
                continue
            rule = str(issue.get("rule_id") or "rule")
            message = str(issue.get("message") or "")
            issue_lines.append(f"- `{rule}`: {message}")
    generated = semantic.get("generated_inputs") or []
    generated_lines = []
    if isinstance(generated, list):
        for item in generated:
            if not isinstance(item, dict):
                continue
            path = str(item.get("path") or "generated input")
            route = str(item.get("route") or "").strip()
            generated_lines.append(f"- `{path}` {route}".rstrip())
    body = [
        "",
        "",
        "runtime semantic gate:",
        f"- verdict: `{verdict}`",
        f"- failed_rule_ids: `{failed_text}`",
    ]
    if issue_lines:
        body.append("- issues:")
        body.extend(issue_lines)
    if generated_lines:
        body.append("- generated input evidence:")
        body.extend(generated_lines[:3])
    return "\n".join(body)


def _format_synthesis_exception(exc: Exception) -> str:
    message = str(exc) or exc.__class__.__name__
    if "MLX runtime requires Apple Silicon/Metal and mlx-lm" in message:
        return (
            "The active local provider is configured for MLX, but this Python "
            "environment cannot load `mlx_lm` or Apple Metal. Start the TUI "
            "with the MLX-enabled interpreter, for example "
            "`/Users/hongjiseung/developer/chemsmart/.venv-mlx/bin/python -m "
            "chemsmart.cli.main agent`, or configure an API frontier provider "
            "in `agent.yaml`. Active-env install command: "
            "`python -m pip install 'mlx-lm==0.31.3'`."
        )
    if "local provider" in message.lower() or "mlx" in message.lower():
        return (
            "The local model provider could not start.\n\n"
            f"{message}\n\n"
            "Use `/doctor` to inspect the active provider, run the TUI from the "
            "MLX-enabled environment, or configure an API frontier provider "
            "in `agent.yaml`."
        )
    return message
