"""Workspace project YAML write approval flow."""

from __future__ import annotations

import shlex
from pathlib import Path

from chemsmart.agent.tui.chat_helpers import (
    _find_project_yaml_candidate_for_write,
)
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    ProjectWriteOverlay,
    ProjectWriteResult,
)
from chemsmart.settings.workspace_project import workspace_project_path


class ProjectWriteMixin:
    """Validate user intent before writing workspace project YAML."""

    def _handle_project_write_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before writing project YAML.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid write-project command", str(exc))
            return

        confirmed = False
        overwrite = False
        project_name: str | None = None
        for part in parts:
            value = part.strip()
            lowered = value.lower()
            if lowered in {"yes", "y", "--yes"}:
                confirmed = True
            elif lowered in {"overwrite", "--overwrite"}:
                overwrite = True
            elif project_name is None:
                project_name = value
            else:
                self.post_error(
                    "Invalid write-project command",
                    "Usage: /write-project [name] (confirmation opens in TUI)",
                )
                return

        candidate = _find_project_yaml_candidate_for_write(
            self.session_root,
            preferred_session_dir=self.current_session_dir(),
        )
        if candidate is None:
            self.post_error(
                "Project write unavailable",
                "Run /project and build a validated project YAML first.",
            )
            return
        if project_name:
            candidate["project_name"] = project_name

        if confirmed:
            candidate["overwrite"] = overwrite
            self._start_project_yaml_write(candidate)
            return
        if self.app.plain:
            self.post_error(
                "Confirmation required",
                "Plain mode uses /write-project [name] yes [overwrite].",
            )
            return
        self._prompt_project_yaml_write(
            candidate,
            prefer_active=project_name is None,
        )

    def _prompt_project_yaml_write(
        self,
        candidate: dict[str, object],
        *,
        prefer_active: bool,
    ) -> None:
        program = str(candidate.get("program") or "gaussian").strip().lower()
        requested = Path(str(candidate.get("project_name") or "project")).stem
        status = self._resolve_workspace_project_status()
        overwrite_project = requested
        if (
            prefer_active
            and status.loaded
            and status.program == program
            and status.project
        ):
            overwrite_project = status.project
        target = workspace_project_path(overwrite_project, program)
        new_project = self._next_available_project_name(requested, program)
        self._pending_project_write_candidate = dict(candidate)
        self.app.push_screen(
            ProjectWriteOverlay(
                target=target,
                overwrite_project=overwrite_project,
                new_project=new_project,
                target_exists=target.exists(),
            ),
            self._handle_project_write_result,
        )

    def _handle_project_write_result(
        self, result: ProjectWriteResult | None
    ) -> None:
        candidate = self._pending_project_write_candidate
        self._pending_project_write_candidate = None
        if result is None or candidate is None:
            self.query_one(FooterWidget).set_hint(
                "Project YAML write cancelled"
            )
            self.focus_composer()
            return
        candidate["project_name"] = result.project_name
        candidate["overwrite"] = result.action == "overwrite"
        self._start_project_yaml_write(candidate)

    def _start_project_yaml_write(self, candidate: dict[str, object]) -> None:
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Writing project YAML…")
        self._current_worker = self.run_project_yaml_write(candidate)

    def _next_available_project_name(
        self, base_name: str, program: str
    ) -> str:
        base = Path(str(base_name or "project")).stem or "project"
        if not workspace_project_path(base, program).exists():
            return base
        suffix = 2
        while workspace_project_path(f"{base}-{suffix}", program).exists():
            suffix += 1
        return f"{base}-{suffix}"
