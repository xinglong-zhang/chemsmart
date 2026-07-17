"""Decision-log streaming and tool-result presentation for the chat screen."""

from __future__ import annotations

from pathlib import Path
from threading import get_ident

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.services.conversation_memory import ConversationMemory
from chemsmart.agent.tui.services.log_tailer import LogTailer


class DecisionLogMixin:
    def _attach_live_tailer(self) -> None:
        session_dir = self.current_session_dir()
        if session_dir is None:
            return
        log_path = session_dir / "decision_log.jsonl"
        if self._tailer_path == log_path:
            return
        if log_path.exists():
            self._attach_tailer(log_path)
            if self._session_poll_timer is not None:
                self._session_poll_timer.stop()
                self._session_poll_timer = None

    def _attach_tailer(self, log_path: Path) -> None:
        self._stop_tailer()
        self._tailer_path = log_path
        self._tailer = LogTailer(log_path, self._on_log_entry)
        self._tailer.start()

    def _stop_tailer(self) -> None:
        if self._tailer is not None:
            self._tailer.stop()
        self._tailer = None
        self._tailer_path = None

    def _on_log_entry(self, entry: dict) -> None:
        app_thread_id = getattr(self.app, "_thread_id", None)
        if app_thread_id == get_ident():
            self._apply_log_entry(entry)
            return
        self.app.call_from_thread(self._apply_log_entry, entry)

    def _current_entity_snapshot(self) -> dict[str, object] | None:
        session_dir = self.current_session_dir()
        session = self.active_agent_session

        if session_dir is not None:
            log_path = session_dir / "decision_log.jsonl"
            if log_path.exists():
                try:
                    memory = ConversationMemory.from_entries(
                        DecisionLog(log_path).read_all()
                    )
                    entities = memory.entities.model_dump(exclude_none=True)
                    return entities or None
                except Exception:
                    pass

        if session is None:
            return None
        entities = session.conversation_history.entities.model_dump(
            exclude_none=True
        )
        return entities or None
