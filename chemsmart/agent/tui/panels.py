"""Chrome and live panels: wordmark, phase chip, jobs, and DAG."""

from __future__ import annotations

import time
from importlib import metadata
from typing import Any, Mapping

from rich.text import Text
from textual.widgets import Static

from chemsmart.agent.voice import human_state

#: node state -> (glyph, style). Minimal color, one rule: success green,
#: failure red, running warning, everything waiting or deferred muted.
_STATE_GLYPHS: dict[str, tuple[str, str]] = {
    "queued": ("[ ]", "dim"),
    "pending": ("[ ]", "dim"),
    "deferred": ("⊘", "dim"),
    "running": ("[•]", "yellow"),
    "engine_complete": ("[»]", "yellow"),
    "validated": ("[✓]", "green"),
    "executed": ("[✓]", "green"),
    "completed": ("[✓]", "green"),
    "failed": ("[✗]", "red"),
    "blocked": ("[✗]", "red"),
    "blocked_unsupported": ("⊘", "dim"),
    "skipped": ("⊘", "dim"),
    "ambiguous": ("[?]", "yellow"),
}


def product_version() -> str:
    try:
        return metadata.version("chemsmart")
    except metadata.PackageNotFoundError:  # pragma: no cover - dev checkout
        return "dev"


#: The logo's vibrational-waveform motif, tiled to width. Exactly two
#: sanctioned placements: under the wordmark, and before "Analysis results".
_SPECTRUM_GLYPHS = "▂▄▆█▆▄▂▁"
SPECTRUM_COLORS = ("#51cf66", "#4dd3e8", "#c2255c", "#ff8fab")


def wordmark() -> Text:
    return Text.assemble(
        ("CHEM", "bold magenta"),
        ("SMART", "bold green"),
        (" Agent", ""),
        (f"  v{product_version()}", "dim"),
    )


def spectrum_rule(width: int, *, plain: bool = False) -> Text:
    """The one signature flourish; a plain dim rule without color."""

    width = max(8, int(width))
    if plain:
        return Text("─" * width, style="dim")
    glyphs = (_SPECTRUM_GLYPHS * (width // len(_SPECTRUM_GLYPHS) + 1))[:width]
    segment = max(1, width // len(SPECTRUM_COLORS))
    text = Text()
    for index, color in enumerate(SPECTRUM_COLORS):
        start = index * segment
        end = width if index == len(SPECTRUM_COLORS) - 1 else start + segment
        if start >= width:
            break
        text.append(glyphs[start:end], style=color)
    return text


def phase_chip(phase: str, hint: str, *, dot_style: str = "bold cyan") -> Text:
    return Text.assemble(
        ("● ", dot_style),
        (phase, "bold"),
        ("  ·  ", "dim"),
        (hint, "dim"),
    )


def state_glyph(state: str) -> tuple[str, str]:
    return _STATE_GLYPHS.get(state, ("[ ]", "dim"))


def dag_rows(nodes: Mapping[str, Mapping[str, Any]]) -> Text:
    """Flat glyph rows over the shared node-state mapping.

    ``nodes`` is ordered: calculation nodes in plan order, then analysis
    nodes. Each value carries label/state/detail/kind.
    """

    text = Text()
    for node_id, node in nodes.items():
        glyph, style = state_glyph(str(node.get("state") or "queued"))
        text.append(f" {glyph} ", style=style)
        text.append(
            node_id, style="bold" if node.get("kind") == "calc" else ""
        )
        label = str(node.get("label") or "")
        if label:
            text.append(f" · {label}", style="dim")
        state = human_state(str(node.get("state") or "queued"))
        text.append(f" · {state}", style=style or "dim")
        detail = str(node.get("detail") or "")
        if detail:
            text.append(f"\n     ↳ {detail}", style="dim")
        text.append("\n")
    if not nodes:
        text.append("no workflow nodes to display", style="dim")
    return text


class DagPanel(Static):
    """The plan with every node's live status; toggled with /dag."""

    def refresh_from(self, nodes: Mapping[str, Mapping[str, Any]]) -> None:
        self.update(dag_rows(nodes))


class JobsPanel(Static):
    """What is really running now: node, command, elapsed wall-clock."""

    def refresh_from(
        self,
        nodes: Mapping[str, Mapping[str, Any]],
        *,
        now: float | None = None,
    ) -> None:
        now = time.monotonic() if now is None else now
        text = Text()
        running = [
            (node_id, node)
            for node_id, node in nodes.items()
            if node.get("state") in {"running", "engine_complete"}
        ]
        queued = [
            node_id
            for node_id, node in nodes.items()
            if node.get("state") in {"queued", "pending"}
            and node.get("kind") == "calc"
        ]
        done = sum(
            1
            for node in nodes.values()
            if node.get("state") in {"validated", "executed", "completed"}
        )
        failed = sum(
            1 for node in nodes.values() if node.get("state") == "failed"
        )
        for node_id, node in running:
            started = node.get("started")
            elapsed = f" · {now - started:.0f} s" if started else ""
            stage = (
                "validating result"
                if node.get("state") == "engine_complete"
                else "engine running"
            )
            text.append("▸ ", style="yellow")
            text.append(f"{node_id}", style="bold")
            text.append(f" · {node.get('label') or ''} · {stage}{elapsed}\n")
            argv = str(node.get("argv") or "")
            if argv:
                text.append(f"  $ {argv}\n", style="dim")
        if not running:
            text.append("no engine running\n", style="dim")
        summary = f"done {done}"
        if failed:
            summary += f" · failed {failed}"
        if queued:
            summary += " · queued: " + ", ".join(queued)
        text.append(summary, style="dim")
        self.update(text)


__all__ = [
    "DagPanel",
    "SPECTRUM_COLORS",
    "JobsPanel",
    "dag_rows",
    "phase_chip",
    "product_version",
    "spectrum_rule",
    "state_glyph",
    "wordmark",
]
