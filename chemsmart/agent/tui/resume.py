"""Discover a workspace's story so a session can be resumed.

Everything here reads durable evidence only: the planning runs' public
transcripts and event streams, the workspace review copies, the per-review
decision logs, and the execution directories. Nothing is written, and no
authority is created -- a pending review found here is re-presented for one
fresh human decision through the ordinary durable approval chain.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path

from .runs import RunSummaryV1, list_runs
from chemsmart.agent.voice import human_state

_PRIVATE = ".chemsmart-agent"


@dataclass(frozen=True)
class SessionStoryV1:
    session_name: str
    task: str
    ended: str
    final_prose: str


@dataclass(frozen=True)
class PendingReviewV1:
    status: str  # "pending" | "approved_unexecuted"
    path: Path
    review_sha256: str
    recovery: str = ""


@dataclass(frozen=True)
class RunRowV1:
    name: str
    kind: str
    workflow_id: str
    state: str
    report_path: Path | None


@dataclass(frozen=True)
class WorkspaceStoryV1:
    story: SessionStoryV1 | None
    runs: tuple[RunRowV1, ...]
    pending: PendingReviewV1 | None
    notes: tuple[str, ...] = field(default_factory=tuple)


def has_history(workspace: Path) -> bool:
    root = Path(workspace) / _PRIVATE
    return any(
        (root / name).is_dir() and any((root / name).iterdir())
        for name in ("runs", "executions", "replays", "reviews")
        if (root / name).is_dir()
    )


def find_latest_session(workspace: Path) -> Path | None:
    runs_root = Path(workspace) / _PRIVATE / "runs"
    if not runs_root.is_dir():
        return None
    candidates = sorted(
        item
        for item in runs_root.iterdir()
        if item.is_dir() and item.name.startswith("live-")
    )
    return candidates[-1] if candidates else None


def _terminal_words(events_file: Path) -> str:
    ended = ""
    reason = ""
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
                    payload = event.get("payload") or {}
                    ended = human_state(str(payload.get("terminal_state")))
                    reason = str(payload.get("reason") or "")
    except OSError:
        return "unreadable"
    if not ended:
        return "ended without a terminal record"
    return f"{ended} ({reason})" if reason and reason != ended else ended


def read_session_story(run_directory: Path) -> SessionStoryV1 | None:
    transcripts = sorted(run_directory.glob("public-transcript-*.json"))
    if not transcripts:
        return None
    try:
        payload = json.loads(transcripts[-1].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    messages = payload.get("transcript") or []
    task = ""
    if len(messages) > 1 and messages[1].get("role") == "user":
        try:
            context = json.loads(messages[1].get("content") or "{}")
            task = str(context.get("task") or "")
        except json.JSONDecodeError:
            task = ""
    final_prose = ""
    for message in reversed(messages):
        if message.get("role") == "assistant" and message.get("content"):
            final_prose = str(message["content"])
            break
    return SessionStoryV1(
        session_name=run_directory.name,
        task=task,
        ended=_terminal_words(run_directory / "events.jsonl"),
        final_prose=final_prose,
    )


def _review_sha256_of(path: Path) -> str:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    record = payload.get("workflow_execution_review", payload)
    return str(record.get("review_sha256") or "")


def find_pending_review(workspace: Path) -> PendingReviewV1 | None:
    from chemsmart.agent.live_session import spent_workflow_approval_ids

    workspace = Path(workspace)
    reviews_root = workspace / _PRIVATE / "reviews"
    if not reviews_root.is_dir():
        return None
    candidates = sorted(
        reviews_root.glob("*.json"),
        key=lambda item: item.stat().st_mtime,
        reverse=True,
    )
    for path in candidates:
        review_sha256 = _review_sha256_of(path)
        if not review_sha256:
            continue
        if spent_workflow_approval_ids(workspace, review_sha256):
            continue  # decided and executed; nothing pending
        scope = workspace / _PRIVATE / "decisions" / review_sha256[:16]
        decisions = scope / "decisions.jsonl"
        approved = False
        if decisions.exists():
            recorded = decisions.read_text(encoding="utf-8")
            approved = (
                '"decision":"approve"' in recorded
                or '"decision": "approve"' in recorded
            )
        if approved:
            if (scope / "bundle.json").exists():
                recovery = (
                    "chemsmart agent run "
                    f"--approval-file {scope / 'bundle.json'} "
                    f"--workspace {workspace} "
                    f"--run-directory {scope / 'run'}"
                )
                return PendingReviewV1(
                    status="approved_unexecuted",
                    path=path,
                    review_sha256=review_sha256,
                    recovery=recovery,
                )
            continue
        return PendingReviewV1(
            status="pending", path=path, review_sha256=review_sha256
        )
    return None


def _workflow_id_of(summary: RunSummaryV1) -> str:
    events = summary.directory / "events.jsonl"
    try:
        with events.open("r", encoding="utf-8") as handle:
            first = handle.readline().strip()
        session_id = str(json.loads(first).get("session_id") or "")
    except (OSError, json.JSONDecodeError):
        return ""
    prefix = "execute-"
    return session_id[len(prefix):] if session_id.startswith(prefix) else ""


def run_rows(workspace: Path) -> tuple[RunRowV1, ...]:
    rows = []
    for summary in list_runs(Path(workspace)):
        state = human_state(summary.terminal_state)
        if summary.report_path is not None:
            state = "completed · report available"
        rows.append(
            RunRowV1(
                name=summary.name,
                kind=summary.kind,
                workflow_id=_workflow_id_of(summary),
                state=state,
                report_path=summary.report_path,
            )
        )
    return tuple(rows)


def workspace_story(workspace: Path) -> WorkspaceStoryV1:
    workspace = Path(workspace)
    latest = find_latest_session(workspace)
    story = read_session_story(latest) if latest is not None else None
    return WorkspaceStoryV1(
        story=story,
        runs=run_rows(workspace),
        pending=find_pending_review(workspace),
    )


__all__ = [
    "PendingReviewV1",
    "RunRowV1",
    "SessionStoryV1",
    "WorkspaceStoryV1",
    "find_latest_session",
    "find_pending_review",
    "has_history",
    "read_session_story",
    "run_rows",
    "workspace_story",
]
