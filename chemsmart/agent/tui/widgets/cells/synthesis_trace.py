"""Public synthesis trace cells for CLI-first ask mode."""

from __future__ import annotations

from rich.console import Group
from rich.markdown import Markdown
from rich.table import Table
from rich.text import Text
from textual import events

from .base import BaseCell


class SynthesisTraceCell(BaseCell):
    """Render observable synthesis steps without exposing hidden reasoning."""

    def __init__(
        self,
        *,
        provider_type: str,
        model: str,
        mode: str,
        status: str,
        command: str = "",
        semantic: dict[str, object] | None = None,
        decision_trace: dict[str, object] | None = None,
        artifact_dir: str = "",
    ) -> None:
        self.source_text = _format_trace_text(
            provider_type=provider_type,
            model=model,
            mode=mode,
            status=status,
            command=command,
            semantic=semantic,
            decision_trace=decision_trace,
            artifact_dir=artifact_dir,
        )
        super().__init__(
            _render_trace(
                provider_type=provider_type,
                model=model,
                mode=mode,
                status=status,
                command=command,
                semantic=semantic,
                decision_trace=decision_trace,
                artifact_dir=artifact_dir,
            ),
            title="Agent Trace",
            classes="agent-cell synthesis-trace-cell",
        )


class FinalAnswerCell(BaseCell):
    """Bottom-of-turn cell containing the final command or answer."""

    def __init__(self, text: str, *, title: str) -> None:
        self.source_text = text
        super().__init__(
            Markdown(text or "_No final output was generated._"),
            title=title,
            classes="agent-cell final-answer-cell",
        )
        self.tooltip = "Click to select and copy this response"

    def on_click(self, event: events.Click) -> None:
        event.stop()
        self.post_message(
            self.CopyRequested(
                self.source_text, self.border_title or "Response"
            )
        )


def _render_trace(
    *,
    provider_type: str,
    model: str,
    mode: str,
    status: str,
    command: str,
    semantic: dict[str, object] | None,
    decision_trace: dict[str, object] | None,
    artifact_dir: str,
) -> Group:
    headline = Text.assemble(
        ("synthesis path: ", "dim"),
        (_lane_label(provider_type), "bold cyan"),
        (" -> ", "dim"),
        ("deterministic CLI harness", "bold"),
        (" -> ", "dim"),
        (status or "unknown", _status_style(status)),
    )
    meta = Text.assemble(
        ("mode ", "dim"),
        (mode or "ask", "cyan"),
        (" · model ", "dim"),
        (model or "auto", "cyan"),
    )
    if artifact_dir:
        meta.append(" · log ", style="dim")
        meta.append(artifact_dir, style="green")

    steps = Table.grid(padding=(0, 2))
    steps.add_column("Step", style="cyan", no_wrap=True)
    steps.add_column("Tool", style="bold", no_wrap=True)
    steps.add_column("Evidence", overflow="fold")

    steps.add_row("1", "yaml_check", "workspace project YAML status checked")
    if _is_api(provider_type):
        action = ""
        if decision_trace:
            action = str(decision_trace.get("action") or "")
        steps.add_row(
            "2",
            "api_intent_router",
            action or "frontier model selected the turn route",
        )
        model_step = "api_command_or_answer"
    else:
        steps.add_row(
            "2",
            "local_command_model",
            "local model used only for CLI command synthesis",
        )
        model_step = "compact_spec_to_command"
    steps.add_row("3", model_step, _command_label(command))
    steps.add_row(
        "4",
        "deterministic_command_parser",
        "real chemsmart CLI command was parsed into user-facing fields",
    )
    steps.add_row("5", "safe_runtime_semantic_gate", _semantic_label(semantic))
    steps.add_row(
        "6", "final_response", "final command or answer rendered last"
    )

    parts: list[object] = [headline, meta, steps]
    public_reasoning = _public_reasoning(decision_trace)
    if public_reasoning:
        reason_table = Table.grid(padding=(0, 1))
        reason_table.add_column(
            "Public reasoning", style="magenta", no_wrap=True
        )
        reason_table.add_column("Evidence", overflow="fold")
        for index, item in enumerate(public_reasoning, start=1):
            reason_table.add_row(f"{index}.", item)
        parts.append(reason_table)
    parts.append(
        Text(
            "Public trace only: hidden model chain-of-thought is not displayed.",
            style="dim",
        )
    )
    return Group(*parts)


def _format_trace_text(
    *,
    provider_type: str,
    model: str,
    mode: str,
    status: str,
    command: str,
    semantic: dict[str, object] | None,
    decision_trace: dict[str, object] | None,
    artifact_dir: str,
) -> str:
    lines = [
        "agent synthesis trace:",
        f"- lane: {_lane_label(provider_type)}",
        f"- mode: {mode or 'ask'}",
        f"- model: {model or 'auto'}",
        f"- status: {status or 'unknown'}",
    ]
    if artifact_dir:
        lines.append(f"- artifact dir: {artifact_dir}")
    if command:
        lines.append(f"- command: `{command}`")
    lines.append(f"- semantic gate: {_semantic_label(semantic)}")
    if decision_trace:
        lines.append(f"- router: {decision_trace.get('router') or 'unknown'}")
        lines.append(
            f"- routed action: {decision_trace.get('action') or 'unknown'}"
        )
        reasoning = _public_reasoning(decision_trace)
        if reasoning:
            lines.append("- public reasoning:")
            lines.extend(f"  - {item}" for item in reasoning)
    lines.append(
        "- note: public trace only; hidden chain-of-thought is not displayed."
    )
    return "\n".join(lines)


def _lane_label(provider_type: str) -> str:
    return "local synthesis" if provider_type == "local" else "api synthesis"


def _is_api(provider_type: str) -> bool:
    return provider_type not in {"", "local", "offline"}


def _command_label(command: str) -> str:
    return command if command else "no chemsmart command emitted"


def _semantic_label(semantic: dict[str, object] | None) -> str:
    if not semantic:
        return "not run for this turn"
    verdict = str(semantic.get("verdict") or "unknown")
    failed = semantic.get("failed_rule_ids") or []
    if isinstance(failed, list) and failed:
        failed_text = ", ".join(str(item) for item in failed)
        return f"{verdict}; failed rules: {failed_text}"
    return f"{verdict}; failed rules: none"


def _status_style(status: str) -> str:
    if status == "ready":
        return "bold green"
    if status in {"needs_clarification", "informational"}:
        return "bold yellow"
    if status in {"infeasible", "reject"}:
        return "bold red"
    return "bold"


def _public_reasoning(trace: dict[str, object] | None) -> list[str]:
    if not isinstance(trace, dict):
        return []
    value = trace.get("reasoning") or trace.get("evidence")
    if not isinstance(value, list):
        return []
    return [str(item) for item in value if str(item).strip()]
