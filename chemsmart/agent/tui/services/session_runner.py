"""Textual worker helpers for running AgentSession methods."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from textual import work

from chemsmart.agent.core import AgentSession


class SessionRunnerMixin:
    session_root: Path
    active_agent_session: AgentSession | None = None
    active_resume_id: str | None = None

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-run",
    )
    def run_agent_session(self, request: str) -> dict[str, Any]:
        self.active_resume_id = None
        self.active_agent_session = AgentSession(
            session_root=str(self.session_root)
        )
        return self.active_agent_session.run(
            request,
            dry_submit=True,
            pause_before_risky=True,
        )

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-build",
    )
    def run_build_session(self, request: str) -> dict[str, Any]:
        """Run a request in project-YAML build mode (tool-loop harness)."""
        self.active_resume_id = None
        self.active_agent_session = AgentSession(
            session_root=str(self.session_root),
            stage_prompt="project_yaml_build.md",
        )
        return self.active_agent_session.run(
            request,
            dry_submit=True,
            pause_before_risky=True,
        )

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-resume",
    )
    def resume_agent_session(
        self,
        session_id: str,
        *,
        cwd_override: str | None = None,
    ) -> dict[str, Any]:
        self.active_resume_id = session_id
        self.active_agent_session = None
        return AgentSession.resume(
            session_id,
            session_root=str(self.session_root),
            dry_submit=True,
            pause_before_risky=True,
            cwd_override=cwd_override,
        )

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-execute",
    )
    def execute_agent_session(self, session_id: str) -> dict[str, Any]:
        self.active_resume_id = session_id
        self.active_agent_session = None
        return AgentSession.resume(
            session_id,
            session_root=str(self.session_root),
            dry_submit=False,
            pause_before_risky=False,
        )

    def current_session_dir(self) -> Path | None:
        if self.active_agent_session is not None:
            return self.active_agent_session.session_dir
        if self.active_resume_id:
            return self.session_root / self.active_resume_id
        return None
