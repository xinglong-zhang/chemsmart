"""Browse the workspace's past executions from their durable evidence."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

_REPORT_RELATIVE = Path("analysis") / "completed-analysis-report.md"


@dataclass(frozen=True)
class RunSummaryV1:
    name: str
    kind: str  # "execution" | "replay"
    directory: Path
    terminal_state: str
    report_path: Path | None


def _terminal_state(events_file: Path) -> str:
    if not events_file.exists():
        return "no event stream"
    state = "in progress"
    try:
        with events_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if event.get("kind") == "runtime_terminated":
                    state = str(
                        (event.get("payload") or {}).get("terminal_state")
                        or "terminated"
                    )
    except OSError:
        return "unreadable"
    return state


def _summarize(directory: Path, kind: str) -> RunSummaryV1:
    # A replay scope keeps its run under <scope>/run; an execution directory
    # is itself the run directory.
    candidates = (directory, directory / "run")
    run_directory = next(
        (
            candidate
            for candidate in candidates
            if (candidate / "events.jsonl").exists()
        ),
        directory,
    )
    report = run_directory / _REPORT_RELATIVE
    return RunSummaryV1(
        name=directory.name,
        kind=kind,
        directory=run_directory,
        terminal_state=_terminal_state(run_directory / "events.jsonl"),
        report_path=report if report.exists() else None,
    )


def list_runs(workspace: Path) -> tuple[RunSummaryV1, ...]:
    """Newest-first summaries of every execution and replay in a workspace."""

    root = Path(workspace) / ".chemsmart-agent"
    found: list[tuple[float, RunSummaryV1]] = []
    for kind, subdir in (("execution", "executions"), ("replay", "replays")):
        base = root / subdir
        if not base.is_dir():
            continue
        for directory in base.iterdir():
            if not directory.is_dir():
                continue
            found.append(
                (directory.stat().st_mtime, _summarize(directory, kind))
            )
    found.sort(key=lambda item: item[0], reverse=True)
    return tuple(summary for _mtime, summary in found)


__all__ = ["RunSummaryV1", "list_runs"]
