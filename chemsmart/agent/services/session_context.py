"""Session ledger, turn state, artifacts, and runtime-controller context."""

from __future__ import annotations

import json
import os
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Protocol

from chemsmart.agent.handles import HandleStore, store_result_handle
from chemsmart.agent.models import CriticVerdict, SessionState, utc_now_iso
from chemsmart.agent.permissions import PermissionPolicy
from chemsmart.agent.runtime.contracts import RuntimeV2Mode
from chemsmart.agent.runtime.orchestrator import RuntimeController
from chemsmart.agent.services.conversation_memory import ConversationMemory
from chemsmart.agent.services.result_codec import json_safe, preview_value
from chemsmart.agent.services.session_store import load_current_session_state

UTC = timezone.utc


class DecisionLog:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(
        self, kind: str, payload: dict[str, Any], rationale: str = ""
    ) -> None:
        entry = {
            "ts": datetime.now(UTC).isoformat(),
            "kind": kind,
            "payload": json_safe(payload),
            "rationale": rationale,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")

    def read_all(self) -> list[dict[str, Any]]:
        if not self.path.exists():
            return []
        with self.path.open(encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]


class SessionContextHost(Protocol):
    session_root: Path
    session_dir: Path | None
    state: SessionState | None
    handle_store: HandleStore | None
    decision_log: DecisionLog | None
    conversation_history: ConversationMemory
    registry: Any
    runtime_v2_mode: RuntimeV2Mode
    _runtime_controller: RuntimeController | None
    _loop_mode_state: tuple[str, bool] | None
    _llm_stats: list[dict[str, Any]]
    _provider: Any | None


class SessionContext:
    def __init__(self, session: SessionContextHost) -> None:
        self.session = session

    def start_new_session(self, request: str) -> None:
        session = self.session
        session_id = new_session_id()
        session.session_dir = session.session_root / session_id
        session.session_dir.mkdir(parents=True, exist_ok=True)
        session.handle_store = HandleStore(session.session_dir)
        session.decision_log = DecisionLog(
            session.session_dir / "decision_log.jsonl"
        )
        session.state = SessionState(
            session_id=session_id,
            cwd=os.path.abspath(os.getcwd()),
            request_started_at=utc_now_iso(),
            turn_index=1,
            current_step_index=0,
            request=request,
            env_snapshot=env_snapshot(),
        )
        session.conversation_history = ConversationMemory()

    def start_new_turn(self, request: str) -> None:
        session = self.session
        assert session.state is not None
        assert session.decision_log is not None
        if self.logged_summary() is None and session.state.plan is not None:
            raise RuntimeError(
                "Cannot start a new request while the current turn is still "
                "open. Resume or finish the existing turn first."
            )
        session.state.turn_index += 1
        session.state.cwd = os.path.abspath(os.getcwd())
        session.state.request_started_at = utc_now_iso()
        session.state.current_step_index = 0
        session.state.total_steps_planned = 0
        session.state.plan = None
        session.state.request = request
        session.state.request_intent = "unknown"
        session.state.env_snapshot = env_snapshot()

    def load_existing(self, session_id: str) -> None:
        session = self.session
        session.session_dir = session.session_root / session_id
        session.state = load_current_session_state(
            session.session_dir, required=True
        )
        session.handle_store = HandleStore(session.session_dir)
        session.decision_log = DecisionLog(
            session.session_dir / "decision_log.jsonl"
        )
        self.refresh_history()
        expected = max(1, len(session.conversation_history.turns))
        if session.state.turn_index != expected:
            session.state.turn_index = expected
            self.save_state()

    def save_state(self) -> None:
        session = self.session
        assert session.session_dir is not None
        assert session.state is not None
        session.state.save(session.session_dir / "session.json")
        session.state.save(session.session_dir / "state.json")

    def refresh_history(self) -> None:
        session = self.session
        if session.decision_log is None:
            session.conversation_history = ConversationMemory()
            return
        session.conversation_history = ConversationMemory.from_entries(
            session.decision_log.read_all()
        )

    def current_turn_entries(self) -> list[dict[str, Any]]:
        session = self.session
        if session.decision_log is None or session.state is None:
            return []
        return session.conversation_history.entries_for_turn(
            session.decision_log.read_all(), session.state.turn_index
        )

    def write_result_artifact(self, step_index: int, result: Any) -> Path:
        session = self.session
        assert session.session_dir is not None
        assert session.state is not None
        path = session.session_dir / (
            f"turn_{session.state.turn_index:02d}_step_{step_index + 1:02d}.json"
        )
        path.write_text(
            json.dumps(json_safe(result), indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return path

    def artifact_paths(self) -> list[Path] | None:
        session = self.session
        assert session.session_dir is not None
        results = [
            entry
            for entry in self.current_turn_entries()
            if entry.get("kind") == "tool_result"
        ]
        paths = []
        for entry in results:
            name = (entry.get("payload") or {}).get("artifact")
            if isinstance(name, str) and name:
                paths.append(session.session_dir / name)
        return paths or None

    def load_completed_results(self) -> list[Any]:
        session = self.session
        assert session.session_dir is not None
        assert session.state is not None
        paths = self.artifact_paths()
        if paths is None:
            paths = [
                session.session_dir / f"step_{index + 1:02d}.json"
                for index in range(session.state.current_step_index)
            ]
        results = []
        for path in paths:
            with path.open(encoding="utf-8") as handle:
                results.append(json.load(handle))
        return results

    def store_result_handle(self, tool_name: str, result: Any) -> str | None:
        return store_result_handle(
            self.session.handle_store,
            tool_name,
            result,
            summary=preview_value(result),
        )

    def tools_called(self) -> list[str]:
        tools: list[str] = []
        for entry in self.current_turn_entries():
            if entry.get("kind") not in {
                "tool_call",
                "tool_preview",
                "tool_use_request",
            }:
                continue
            tool = (entry.get("payload") or {}).get("tool")
            if isinstance(tool, str) and tool not in tools:
                tools.append(tool)
        return tools

    def collect_prior_results(
        self, results: list[Any], tool_name: str
    ) -> list[Any]:
        state = self.session.state
        assert state is not None and state.plan is not None
        return [
            result
            for step, result in zip(state.plan.steps, results)
            if step.tool == tool_name
        ]

    def find_prior_result(self, results: list[Any], tool_name: str) -> Any:
        matches = self.collect_prior_results(results, tool_name)
        return matches[0] if matches else None

    def logged_verdict(self) -> CriticVerdict | None:
        entries = [
            entry
            for entry in self.current_turn_entries()
            if entry.get("kind") == "critic_verdict"
        ]
        return (
            CriticVerdict.model_validate(entries[-1]["payload"])
            if entries
            else None
        )

    def logged_summary(self) -> dict[str, Any] | None:
        entries = [
            entry
            for entry in self.current_turn_entries()
            if entry.get("kind") == "session_summary"
        ]
        payload = entries[-1].get("payload") if entries else None
        return payload if isinstance(payload, dict) else None

    def prompt_meta(self, **extra: Any) -> dict[str, Any]:
        session = self.session
        meta: dict[str, Any] = dict(extra)
        if session.state is None:
            return meta
        meta.update(
            {
                "session_id": session.state.session_id,
                "turn_index": session.state.turn_index,
                "request_intent": session.state.request_intent,
            }
        )
        if session._loop_mode_state is not None:
            meta["approval_mode"], meta["yolo"] = session._loop_mode_state
        controller = session._runtime_controller
        if controller is not None:
            if controller.state.active_project is not None:
                meta["active_project"] = controller.state.active_project.name
            if controller.state.previous_command:
                meta["previous_command"] = controller.state.previous_command
        return meta

    def ensure_runtime_controller(self) -> RuntimeController | None:
        session = self.session
        if session.runtime_v2_mode is RuntimeV2Mode.OFF:
            return None
        assert session.state is not None
        assert session.session_dir is not None
        if session._runtime_controller is None:
            session._runtime_controller = RuntimeController(
                session_dir=session.session_dir,
                session_id=session.state.session_id,
                registry=session.registry,
                mode=session.runtime_v2_mode,
            )
        return session._runtime_controller

    def runtime_metadata(self) -> dict[str, Any]:
        session = self.session
        controller = session._runtime_controller
        if controller is None:
            return {"mode": session.runtime_v2_mode.value}
        selection = controller.selection
        return {
            "mode": session.runtime_v2_mode.value,
            "phase": controller.state.phase.value,
            "exposed_tools": list(selection.direct) if selection else [],
            "shadow_violations": list(controller.state.shadow_violations),
            "event_log": str(controller.store.path),
            "state_snapshot": str(
                controller.session_dir / "runtime_state.json"
            ),
        }

    def has_pending_ask_user(self) -> bool:
        state = self.session.state
        if state is None or self.logged_summary() is not None:
            return False
        return bool(state.pending_ask_user and state.pending_messages)

    def resume_pending_ask_user(self, answer: str) -> None:
        session = self.session
        assert session.state is not None
        assert session.decision_log is not None
        session.state.cwd = os.path.abspath(os.getcwd())
        session.state.env_snapshot = env_snapshot()
        pending = session.state.pending_ask_user or {}
        session.decision_log.write(
            "ask_user_answer",
            {
                "question": pending.get("question"),
                "options": pending.get("options") or [],
                "answer": answer,
            },
            rationale=answer,
        )

    def log_loop_mode(self, policy: PermissionPolicy) -> None:
        session = self.session
        assert session.decision_log is not None
        previous = session._loop_mode_state
        session.decision_log.write(
            "mode_change",
            {
                "from_mode": previous[0] if previous else None,
                "to_mode": policy.mode.value,
                "yolo": policy.yolo,
            },
        )
        session._loop_mode_state = (policy.mode.value, policy.yolo)

    def provider_name(self) -> str:
        session = self.session
        if session._llm_stats:
            return session._llm_stats[-1].get("provider_name") or "unknown"
        return getattr(session._provider, "name", None) or "unknown"

    def resolved_model(self) -> str:
        session = self.session
        if session._llm_stats:
            return session._llm_stats[-1].get("resolved_model") or "unknown"
        return getattr(session._provider, "default_model", None) or "unknown"

    def total_tokens(self, field: str) -> int:
        return sum(
            int(stat.get(field) or 0) for stat in self.session._llm_stats
        )


def new_session_id() -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def env_snapshot() -> dict[str, str | None]:
    return {
        "AI_PROVIDER": os.environ.get("AI_PROVIDER"),
        "PWD": os.path.abspath(os.getcwd()),
    }


__all__ = ["DecisionLog", "SessionContext", "env_snapshot", "new_session_id"]
