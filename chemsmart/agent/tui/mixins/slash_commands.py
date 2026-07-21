"""Slash command routing and command-local presenters."""

from __future__ import annotations

import threading
from typing import Callable, Iterable

from click.testing import CliRunner
from rich.table import Table

from chemsmart.agent.cli_commands import agent
from chemsmart.agent.permissions import ApprovalDecision
from chemsmart.agent.services.cli_presenters import sanitize_inline_cli_output
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.screens.sessions import SessionsScreen
from chemsmart.agent.tui.services.session_index import agent_session_dirs
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CriticVerdictCell,
    PlanCell,
)
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import PermissionModeOverlay
from chemsmart.agent.tui.widgets.transcript import Transcript


class SlashCommandsMixin:
    """Dispatch deterministic slash commands to focused handlers."""

    def _handle_slash_command(self, raw: str) -> None:
        command, _, remainder = raw.partition(" ")
        argument = remainder.strip()
        handler = self._slash_command_handlers().get(command)
        if handler is None:
            self.post_error("Unknown command", raw)
            return
        handler(argument)

    def _slash_command_handlers(self) -> dict[str, Callable[[str], None]]:
        return {
            "/help": lambda _arg: self._show_help(),
            "/quit": self._slash_quit,
            "/exit": self._slash_quit,
            "/clear": self._slash_clear,
            "/sessions": self._slash_sessions,
            "/resume": self._slash_resume,
            "/jobs": lambda _arg: self.action_show_calculations(),
            "/runs": lambda _arg: self.action_show_calculations(),
            "/queue": lambda _arg: self._show_queue_snapshot(),
            "/server": self._switch_active_server,
            "/cancel": self._slash_cancel,
            "/extract": self._slash_extract,
            "/molecule": self._slash_molecule,
            "/tools": lambda _arg: self._run_inline_cli(
                ["tools"], title="Tools"
            ),
            "/doctor": lambda _arg: self._run_inline_cli(
                ["doctor"], title="Doctor"
            ),
            "/init": self._handle_init_command,
            "/write-project": self._handle_project_write_command,
            "/wizard": self._handle_wizard_probe_command,
            "/wizard-verify": self._handle_wizard_verify_command,
            "/wizard-refresh": self._handle_wizard_refresh_command,
            "/wizard-write": self._handle_wizard_write_command,
            "/dryrun": self._slash_dryrun,
            "/critic": self._slash_critic,
            "/plan": self._slash_plan,
            "/rationale": self._slash_rationale,
            "/allow": lambda _arg: self._resolve_pending_approval(
                ApprovalDecision.ALLOW_ONCE
            ),
            "/allow-session": lambda _arg: self._resolve_pending_approval(
                ApprovalDecision.ALLOW_SESSION
            ),
            "/deny": lambda _arg: self._resolve_pending_approval(
                ApprovalDecision.DENY
            ),
            "/permissions": self._slash_permissions,
            "/yolo": self._slash_yolo,
            "/execute": self._handle_execute_command,
            "/submit": lambda arg: self._handle_ready_command_execution(
                "sub", arg
            ),
            "/run": lambda arg: self._handle_ready_command_execution(
                "run", arg
            ),
        }

    def _slash_quit(self, _argument: str) -> None:
        self._reset_request_state(clear_transcript=False, clear_session=True)
        self.app.exit()

    def _slash_clear(self, _argument: str) -> None:
        if not self._guard_phase("/clear", {Phase.IDLE, Phase.FINISHED}):
            return
        self._stop_tailer()
        self._reset_request_state(clear_transcript=True, clear_session=True)
        self.notify("Transcript cleared.", timeout=3)
        self.focus_composer()

    def _slash_sessions(self, _argument: str) -> None:
        if not self._guard_phase("/sessions", {Phase.IDLE}):
            return
        if self.app.plain:
            self._show_sessions_snapshot()
            return
        self.app.push_screen(SessionsScreen(self.session_root))

    def _slash_resume(self, argument: str) -> None:
        if not self._guard_phase("/resume", {Phase.IDLE}):
            return
        if argument:
            self._resume_or_prompt(argument)
            return
        self._slash_sessions("")

    def _slash_cancel(self, argument: str) -> None:
        job_id, confirmed = self._parse_cancel_argument(argument)
        if not job_id:
            self.post_error("Missing job id", "Usage: /cancel <job-id> [yes]")
            return
        if confirmed:
            self._handle_cancel_confirmation(job_id, "yes")
            return
        if self.app.plain:
            self.post_error(
                "Confirmation required",
                "Plain mode uses /cancel <job-id> yes.",
            )
            return
        self._confirm_cancel(job_id)

    def _slash_extract(self, argument: str) -> None:
        if argument:
            self._extract_job_result(argument)
            return
        self.post_error("Missing target", "Usage: /extract <job-id|inputfile>")

    def _slash_molecule(self, argument: str) -> None:
        if argument:
            self._show_molecule(argument)
            return
        self.post_error("Missing path", "Usage: /molecule <path>")

    def _slash_dryrun(self, _argument: str) -> None:
        if not self._guard_phase("/dryrun", {Phase.DRY_RUN_READY}):
            return
        if self._current_request:
            self.start_request(self._current_request)
            return
        self.post_error("No request", "There is no active request.")

    def _slash_critic(self, _argument: str) -> None:
        if not self._guard_phase(
            "/critic", {Phase.PLANNING, Phase.DRY_RUN_READY}
        ):
            return
        if self._current_verdict is None:
            self.post_error(
                "No critic verdict", "The critic has not finished yet."
            )
            return
        self.query_one(Transcript).add_cell(
            CriticVerdictCell(self._current_verdict)
        )

    def _slash_plan(self, _argument: str) -> None:
        allowed = {Phase.PLANNING, Phase.DRY_RUN_READY, Phase.RUNNING}
        if not self._guard_phase("/plan", allowed):
            return
        if not self._current_plan_text:
            self.post_error("No plan", "The planner has not finished yet.")
            return
        self.query_one(Transcript).add_cell(PlanCell(self._current_plan_text))

    def _slash_rationale(self, _argument: str) -> None:
        rationale = (
            self._current_plan.rationale if self._current_plan else None
        )
        if rationale:
            self.post_agent_message(rationale, title="Rationale")
            return
        self.post_error("No rationale", "No planner rationale is available.")

    def _slash_permissions(self, argument: str) -> None:
        if argument:
            if not self._apply_permission_command(argument):
                self.post_error(
                    "Unknown permissions command",
                    "Use /permissions permission|driving.",
                )
            return
        if self.app.plain:
            self.post_agent_message(
                "```\n"
                f"mode: {self._permission_mode.value}\n"
                f"yolo: {'on' if self._yolo_enabled else 'off'}\n"
                "Use /permissions permission|driving and /yolo on|off.\n```",
                title="Permissions",
            )
            return
        self.app.push_screen(
            PermissionModeOverlay(
                mode=self._permission_mode, yolo=self._yolo_enabled
            ),
            self._handle_permission_mode_result,
        )

    def _slash_yolo(self, argument: str) -> None:
        if argument:
            self._set_yolo(argument)
            return
        self.post_error("Missing value", "Usage: /yolo on|off")

    def _show_help(self) -> None:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Phase", style="dim", no_wrap=True)
        table.add_column("Command", style="bold")
        table.add_column("Description")
        rows = [
            ("[A]", "/help", "show this help"),
            ("[A]", "/jobs", "open the jobs panel"),
            ("[A]", "/runs · Ctrl+B", "monitor calculations and logs"),
            ("[A]", "/queue", "show the current queue snapshot"),
            ("[A]", "/server <name>", "switch the active HPC server"),
            ("[A]", "/molecule <path>", "load and preview a molecule"),
            ("[A]", "/cancel <job-id> [yes]", "cancel a queued/running job"),
            ("[F]", "/extract <job-id|inputfile>", "parse a final result"),
            ("[D]", "/dryrun", "regenerate the current dry-run"),
            ("[A]", "/allow", "allow the focused tool request once"),
            (
                "[A]",
                "/allow-session",
                "allow the focused tool for the rest of this session",
            ),
            ("[A]", "/deny", "deny the focused tool request"),
            ("[A]", "/permissions", "toggle permission/driving mode"),
            ("[A]", "/yolo on|off", "toggle risky-tool autonomy"),
            ("[F]", "/execute [yes]", "submit the dry-run to HPC for real"),
            ("[F]", "/submit [yes]", "submit validated chemsmart sub command"),
            ("[F]", "/run [yes]", "execute validated chemsmart run command"),
            ("[P,D]", "/critic", "show the current critic verdict"),
            ("[P,D,R]", "/plan", "show the current plan"),
            ("[A]", "/rationale", "show planner rationale"),
            ("[I,F]", "/clear", "clear the transcript"),
            ("[I]", "/sessions", "browse recent sessions"),
            ("[I]", "/resume <session-id>", "load or continue a session"),
            ("[A]", "/tools", "list registered tools"),
            ("[A]", "/wizard-refresh <name> [--force]", "refresh node cache"),
            ("[A]", "/wizard-verify <name>", "verify server transport wiring"),
            ("[A]", "/doctor", "run inline diagnostics"),
            (
                "[F]",
                "/write-project [name]",
                "write or replace the latest validated project YAML",
            ),
            ("[A]", "/quit", "exit the TUI"),
        ]
        for row in rows:
            table.add_row(*row)
        self.query_one(Transcript).add_cell(
            AgentMessageCell(table, title="Help")
        )

    def _handle_init_command(self, argument: str) -> None:
        if (
            self._active_provider_config is not None
            and self._active_provider_config.type == "local"
        ):
            self.post_error(
                "Init unavailable for local provider",
                "Project YAML build mode uses the tool-loop harness, which "
                "needs a tool-calling provider (anthropic/openai).",
            )
            return
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting "
                "project YAML build mode.",
            )
            return
        argument = argument.strip()
        self._build_mode = False
        self._sync_footer_provider()
        if argument:
            self.start_request(
                "Build a workspace project YAML from this reported method. "
                "Validate and critique it before writing; ask for approval "
                f"before write_project_yaml.\n\n{argument}"
            )
            return
        self.post_agent_message(
            (
                "**Project YAML request helper**\n\n"
                "Paste your reported computational method (a paper's "
                "*Computational Details*, a method sentence, or a few facts) "
                "as a normal message, and the unified agent will use the "
                "project-YAML harness to "
                "`extract → render → validate → critique` a chemsmart project "
                "YAML, then `write` it into this workspace's "
                "`.chemsmart/<program>/` folder once you approve.\n\n"
                'Tip: name the project (e.g. "call it co2") and say the '
                "program (Gaussian or ORCA)."
            ),
            title="Init: Project YAML",
        )
        self.query_one(FooterWidget).set_hint("Paste a project-YAML request")
        self.focus_composer()

    def _show_sessions_snapshot(self) -> None:
        lines = ["Sessions", "", "Use /resume <session-id> to load one."]
        if self.session_root.exists():
            session_dirs = agent_session_dirs(self.session_root)
        else:
            session_dirs = []
        if not session_dirs:
            lines.extend(["", "No sessions found."])
        for session_dir in session_dirs[:10]:
            lines.append(f"- {session_dir.name}")
        body = "\n".join(lines)
        self.post_agent_message(f"```\n{body}\n```", title="Sessions")

    def _parse_cancel_argument(self, argument: str) -> tuple[str | None, bool]:
        if not argument:
            return None, False
        parts = argument.split(maxsplit=1)
        if len(parts) == 1:
            return parts[0], False
        if parts[1].strip().lower() in {"y", "yes"}:
            return parts[0], True
        self.post_error(
            "Unknown confirmation",
            "Use /cancel <job-id> yes to confirm.",
        )
        return None, False

    def _handle_inline_approval(self, argument: str) -> bool:
        if not argument:
            return False
        keyword, _, remainder = argument.partition(" ")
        keyword = keyword.strip().lower()
        corrective_text = remainder.strip() or None
        if keyword in {"y", "yes"}:
            self._resolve_pending_approval(ApprovalDecision.ALLOW_ONCE)
            return True
        if keyword in {"n", "no"}:
            self._resolve_pending_approval(ApprovalDecision.DENY)
            return True
        if keyword in {"s", "session"}:
            self._resolve_pending_approval(ApprovalDecision.ALLOW_SESSION)
            return True
        if keyword in {"r", "revise"}:
            if not corrective_text:
                self.post_error(
                    "Missing instruction",
                    "Usage: /run revise <instruction> or /submit revise <instruction>.",
                )
                return True
            self.start_request(self._corrected_request(corrective_text))
            return True
        self.post_error(
            "Unknown approval response",
            "Use yes, no, session, or revise <instruction>.",
        )
        return True

    def _run_inline_cli(self, args: Iterable[str], *, title: str) -> None:

        command_args = list(args)

        def run_inline() -> None:
            try:
                result = CliRunner().invoke(
                    agent, command_args, catch_exceptions=False
                )
            except Exception as exc:
                self.app.call_from_thread(self.post_error, title, str(exc))
                return

            text = sanitize_inline_cli_output(result.output)
            text = text or f"{title} completed."

            def publish() -> None:
                if result.exit_code == 0:
                    self.post_agent_message(f"```\n{text}\n```", title=title)
                else:
                    self.post_error(title, text)

            self.app.call_from_thread(publish)

        threading.Thread(target=run_inline, daemon=True).start()
