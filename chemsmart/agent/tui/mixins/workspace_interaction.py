"""Keyboard, palette, and workspace interactions for the chat screen."""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    select_workspace_project,
)
from chemsmart.agent.synthesis import resolve_default_project
from chemsmart.agent.tui.chat_helpers import _yaml_footer_label
from chemsmart.agent.tui.slash_catalog import (
    SLASH_PALETTE_COMMANDS as _SLASH_PALETTE_COMMANDS,
)
from chemsmart.agent.tui.widgets.cells import AgentMessageCell, ErrorCell
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    HistorySearchOverlay,
    ProjectYamlOverlay,
    ShortcutOverlay,
    ToolActivityOverlay,
)
from chemsmart.agent.tui.widgets.slash_palette import (
    SlashCommandPalette,
    SlashPaletteItem,
)
from chemsmart.agent.tui.widgets.transcript import Transcript
from chemsmart.settings.workspace_project import (
    WorkspaceProjectStatus,
    resolve_workspace_project,
)


class WorkspaceInteractionMixin:
    """Present navigation controls and select workspace project state."""

    def action_show_project_yaml(self) -> None:
        status = self._resolve_workspace_project_status()
        self._workspace_project_status = status
        self.query_one(FooterWidget).set_yaml_status(
            loaded=status.loaded,
            label=_yaml_footer_label(status),
        )
        if not status.candidates:
            if self.app.plain:
                self.post_agent_message(
                    "YAML MISSING\n\nBuild one with /init, then save it with "
                    "/write-project.",
                    title="Workspace YAML",
                )
                return
            self.app.push_screen(
                ProjectYamlOverlay(
                    title="Workspace YAML",
                    candidates=(),
                )
            )
            return
        if self.app.plain:
            if not status.loaded or status.path is None:
                candidates = "\n".join(
                    f"- `{path.parent.name}:{path.stem}`: `{path}`"
                    for path in status.candidates
                )
                self.post_agent_message(
                    "Multiple workspace YAML files are available. Select one "
                    "with an explicit project name in the request or CLI.\n\n"
                    + candidates,
                    title="Workspace YAML selection",
                )
                return
            try:
                yaml_text = status.path.read_text(encoding="utf-8")
            except OSError as exc:
                yaml_text = f"# Failed to read {status.path}: {exc}\n"
            self.post_agent_message(
                f"path: `{status.path}`\n\n```yaml\n{yaml_text.rstrip()}\n```",
                title=f"Workspace YAML: {status.program}:{status.project}",
            )
            return
        self.app.push_screen(
            ProjectYamlOverlay(
                title="Workspace project YAML",
                candidates=status.candidates,
                selected_path=status.path,
            ),
            self._on_project_yaml_selected,
        )

    def _on_project_yaml_selected(self, path: Path | None) -> None:
        if path is not None:
            self._activate_workspace_project_path(path)
        self.focus_composer()

    def action_show_shortcuts(self) -> None:
        if self.app.plain:
            self._show_help()
            return
        config = getattr(self.app, "tui_config", None)
        bindings = config.keybindings if config is not None else {}
        self.app.push_screen(ShortcutOverlay(bindings))

    def action_toggle_transcript(self) -> None:
        transcript = self.query_one(Transcript)
        expanded = transcript.toggle_detail_mode()
        self.query_one(FooterWidget).set_hint(
            "Transcript detail expanded" if expanded else "Transcript compact"
        )

    def action_show_activity(self) -> None:
        entries = [
            self._tool_cells[call_id].activity_snapshot()
            for call_id in self._tool_order
            if call_id in self._tool_cells
        ]
        context = {
            "phase": self.tui_state.phase.value,
            "operation": self.tui_state.operation,
            "progress": self.tui_state.tool_progress or "no active tool",
        }
        if self.app.plain:
            if not entries:
                self.post_agent_message(
                    "No tool activity in this turn.", title="Tool activity"
                )
                return
            lines = [
                f"phase: {context['phase']}",
                f"operation: {context['operation']}",
                f"workflow: {context['progress']}",
                "",
            ]
            lines.extend(
                f"- {entry['status']} {entry['tool']} {entry['elapsed']}"
                for entry in entries
            )
            self.post_agent_message("\n".join(lines), title="Tool activity")
            return
        self.app.push_screen(ToolActivityOverlay(entries, context=context))

    def action_search_history(self) -> None:
        if not self._request_history:
            self.notify("No request history yet.", timeout=2)
            return
        if not self.app.plain:
            self.app.push_screen(
                HistorySearchOverlay(self._request_history),
                self._handle_history_search,
            )
            return
        self._history_cursor = max(0, self._history_cursor - 1)
        composer = self.query_one(Composer)
        composer.load_text(self._request_history[self._history_cursor])
        composer.focus()
        self.query_one(FooterWidget).set_hint(
            f"History {self._history_cursor + 1}/{len(self._request_history)}"
        )

    def _handle_history_search(self, value: str | None) -> None:
        if not value:
            self.focus_composer()
            return
        composer = self.query_one(Composer)
        composer.load_text(value)
        composer.focus()
        self.query_one(FooterWidget).set_hint("History request restored")

    def action_context_tab(self) -> None:
        palette = self.query_one(SlashCommandPalette)
        composer = self.query_one(Composer)
        if palette.is_open:
            selected = palette.selected_item()
            if selected is not None:
                composer.load_text(f"{selected.command} ")
                palette.hide()
            return
        if not self._worker_is_busy():
            return
        draft = composer.resolve_text().strip()
        footer = self.query_one(FooterWidget)
        if draft and self._queued_prompt is None:
            self._queued_prompt = draft
            composer.clear_text()
            footer.update_draft("")
            footer.set_queued_prompt(True)
            footer.set_hint(
                "Follow-up queued · Tab on an empty draft restores it"
            )
            return
        if not draft and self._queued_prompt is not None:
            composer.load_text(self._queued_prompt)
            self._queued_prompt = None
            footer.set_queued_prompt(False)
            footer.set_hint("Queued follow-up restored for editing")
            return
        if draft and self._queued_prompt is not None:
            self.notify(
                "One follow-up is already queued.",
                severity="warning",
                timeout=3,
            )

    def action_soft_cancel(self) -> None:
        if self._quit_armed:
            if self._quit_timer is not None:
                self._quit_timer.stop()
                self._quit_timer = None
            self.app.exit()
            return
        self._quit_armed = True
        self.query_one(FooterWidget).set_hint(
            "Press Ctrl+C again within 3s to quit"
        )
        if self._current_worker and not self._current_worker.is_finished:
            self.post_agent_message(
                "Soft cancel requested. The current step may still finish "
                "because the underlying agent run is not cooperatively "
                "interruptible yet."
            )
        if self._quit_timer is not None:
            self._quit_timer.stop()
        self._quit_timer = self.set_timer(3, self._disarm_soft_cancel)

    def action_quit_if_empty(self) -> None:
        composer = self.query_one(Composer)
        if composer.text.strip():
            return
        self._reset_request_state(clear_transcript=False, clear_session=True)
        self.app.exit()

    def action_refresh_screen(self) -> None:
        self.refresh(repaint=True, layout=True)

    def action_dismiss_overlay(self) -> None:
        if self.app.screen is not self:
            self.app.pop_screen()
            self.call_after_refresh(self.focus_composer)
            return
        self.query_one(SlashCommandPalette).hide()
        self.focus_composer()

    def focus_composer(self) -> None:
        self.query_one(Composer).focus()

    def post_agent_message(self, text: str, *, title: str = "Agent") -> None:
        self.query_one(Transcript).add_cell(
            AgentMessageCell(text, title=title)
        )

    def post_error(
        self, title: str, message: str, details: dict | None = None
    ) -> None:
        self.query_one(Transcript).add_cell(ErrorCell(title, message, details))

    def _update_slash_palette(self, text: str) -> None:
        palette = self.query_one(SlashCommandPalette)
        raw = text.lstrip()
        if not raw.startswith("/") or " " in raw:
            palette.hide()
            return
        query = raw[1:].lower()
        matches = self._slash_palette_items(query)
        palette.show_matches(query=query, matches=matches)

    def _slash_palette_items(self, query: str) -> list[SlashPaletteItem]:
        shortcuts = {
            "/help": (
                getattr(self.app, "tui_config", None).keybindings.get(
                    "show_shortcuts", "f1"
                )
                if getattr(self.app, "tui_config", None) is not None
                else "f1"
            ),
            "/jobs": "ctrl+b",
            "/runs": "ctrl+b",
        }
        items: list[SlashPaletteItem] = []
        for command, description in _SLASH_PALETTE_COMMANDS:
            aliases = {
                "/permissions": ("perms",),
                "/quit": ("q",),
                "/resume": ("continue",),
            }.get(command, ())
            if not (
                command[1:].startswith(query)
                or any(alias.startswith(query) for alias in aliases)
            ):
                continue
            reason = self._slash_unavailable_reason(command)
            items.append(
                SlashPaletteItem(
                    command,
                    description,
                    enabled=reason is None,
                    unavailable_reason=reason or "",
                    shortcut=shortcuts.get(command, ""),
                    aliases=aliases,
                )
            )
        return items

    def _slash_unavailable_reason(self, command: str) -> str | None:
        if (
            command in {"/allow", "/allow-session", "/deny"}
            and not self._pending_approval
        ):
            return "no approval is pending"
        if command == "/execute" and not self._last_dry_run_session_id:
            return "no validated dry-run is ready"
        if command == "/run" and not self._can_approve_or_execute("run"):
            return "no validated local command is ready"
        if command == "/submit" and not self._can_approve_or_execute("sub"):
            return "no validated submission command is ready"
        if command == "/critic" and self._current_verdict is None:
            return "no critic verdict yet"
        if command in {"/plan", "/rationale"} and self._current_plan is None:
            return "no active plan yet"
        return None

    def accept_slash_palette(self) -> str | bool | None:
        palette = self.query_one(SlashCommandPalette)
        if not palette.is_open:
            return None
        selected = palette.selected_item()
        if selected is None:
            reason = next(
                (
                    item.unavailable_reason
                    for item in palette.items
                    if item.unavailable_reason
                ),
                "No available command is selected.",
            )
            self.notify(reason, severity="warning", timeout=3)
            return False
        palette.hide()
        return selected.command

    def _sync_footer_provider(self) -> None:
        config = self._active_provider_config
        self._workspace_project_status = (
            self._resolve_workspace_project_status()
        )
        if config is not None:
            role = "synthesis" if config.type == "local" else "unified"
            self.query_one(FooterWidget).set_provider_model(
                f"{role}:{config.type}",
                config.model,
                project=(
                    self._workspace_project_status.project
                    or config.project
                    or resolve_default_project()
                ),
            )
        self.query_one(FooterWidget).set_yaml_status(
            loaded=self._workspace_project_status.loaded,
            label=_yaml_footer_label(self._workspace_project_status),
        )

    def _resolve_workspace_project_status(self) -> WorkspaceProjectStatus:
        selected_path = self._selected_workspace_project_path
        state_project = current_workflow_state().project
        if selected_path is None and state_project is not None:
            selected_path = Path(state_project.path).resolve()

        status = resolve_workspace_project(selected_path=selected_path)
        if not status.loaded and status.candidates:
            configured = (
                str(self._active_provider_config.project or "").strip()
                if self._active_provider_config is not None
                else ""
            )
            matches = tuple(
                path for path in status.candidates if path.stem == configured
            )
            if len(matches) == 1:
                status = resolve_workspace_project(selected_path=matches[0])

        self._selected_workspace_project_path = (
            status.path if status.loaded else None
        )
        return status

    def _activate_workspace_project_path(self, path: Path) -> bool:
        resolved = path.expanduser().resolve()
        status = resolve_workspace_project(selected_path=resolved)
        if not status.loaded or status.path is None:
            self.post_error(
                "Project selection failed",
                f"{resolved} is not a workspace project YAML.",
            )
            return False
        state_delta = select_workspace_project(
            status.project,
            status.program,
        )
        if not state_delta.get("selected"):
            self.post_error(
                "Project selection failed",
                str(state_delta.get("rule_id") or resolved),
            )
            return False
        self._selected_workspace_project_path = status.path
        self._workspace_project_status = status
        if self.active_synthesis_session is not None:
            self.active_synthesis_session.default_project = status.project
        footer = self.query_one(FooterWidget)
        footer.set_yaml_status(
            loaded=True,
            label=_yaml_footer_label(status),
        )
        footer.set_hint(f"Active project: {status.program}:{status.project}")
        return True

    def _apply_workspace_project_tool_result(
        self,
        tool_name: str,
        result: object,
    ) -> None:
        if tool_name not in {
            "read_project_yaml",
            "write_project_yaml",
            "update_project_yaml",
        } or not isinstance(result, dict):
            return
        raw_path = result.get("written_path") or result.get("path")
        if raw_path:
            self._activate_workspace_project_path(Path(str(raw_path)))
