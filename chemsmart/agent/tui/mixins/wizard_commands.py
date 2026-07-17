"""Wizard slash commands and direct tool lifecycle presentation."""

from __future__ import annotations

import json
import shlex

from rich.table import Table

from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.footer import FooterWidget


class WizardCommandsMixin:
    def _handle_wizard_probe_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard command", str(exc))
            return
        if not parts or len(parts) > 2:
            self.post_error(
                "Invalid wizard command",
                "Usage: /wizard <name> [host]",
            )
            return
        server_name = parts[0]
        ssh_host_hint = parts[1] if len(parts) == 2 else None
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard is probing…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_probe",
            {
                "server_name": server_name,
                "ssh_host_hint": ssh_host_hint,
            },
        )

    def _handle_wizard_verify_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard-verify command", str(exc))
            return
        if len(parts) != 1:
            self.post_error(
                "Invalid wizard-verify command",
                "Usage: /wizard-verify <name>",
            )
            return
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard verify is checking transport wiring…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_verify",
            {"server_name": parts[0]},
        )

    def _handle_wizard_refresh_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard-refresh command", str(exc))
            return
        if not parts or len(parts) > 2:
            self.post_error(
                "Invalid wizard-refresh command",
                "Usage: /wizard-refresh <name> [--force]",
            )
            return
        force = False
        if len(parts) == 2:
            if parts[1] not in {"force", "--force"}:
                self.post_error(
                    "Invalid wizard-refresh command",
                    "Usage: /wizard-refresh <name> [--force]",
                )
                return
            force = True
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard refresh is probing cache state…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_refresh",
            {"server_name": parts[0], "force": force},
        )

    def _handle_wizard_write_command(self, argument: str) -> None:
        if self._current_worker and not self._current_worker.is_finished:
            self.post_error(
                "Session already running",
                "Wait for the current request to finish before starting a new request.",
            )
            return
        try:
            parts = shlex.split(argument)
        except ValueError as exc:
            self.post_error("Invalid wizard-write command", str(exc))
            return
        overwrite = False
        if parts:
            if len(parts) == 1 and parts[0] in {"overwrite", "--overwrite"}:
                overwrite = True
            else:
                self.post_error(
                    "Invalid wizard-write command",
                    "Usage: /wizard-write [overwrite]",
                )
                return
        probe = self._latest_wizard_probe
        if not isinstance(probe, dict):
            self.post_error(
                "Wizard write unavailable",
                "Run /wizard <name> [host] first.",
            )
            return
        validation = probe.get("validation")
        if isinstance(validation, dict) and not validation.get("ok", False):
            self.post_error(
                "Wizard write unavailable",
                "The latest wizard result did not validate.",
            )
            return
        server_name = str(probe.get("server_name") or "")
        yaml_text = str(probe.get("yaml_text") or "")
        if not server_name or not yaml_text:
            self.post_error(
                "Wizard write unavailable",
                "The latest wizard result is incomplete.",
            )
            return
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.PLANNING)
        footer.set_hint("Wizard is writing…")
        self._current_worker = self.run_slash_tool_request(
            "wizard_write",
            {
                "server_name": server_name,
                "yaml_text": yaml_text,
                "overwrite": overwrite,
            },
        )

    def _set_pending_tool_context(
        self,
        description: str,
        arguments: dict,
    ) -> None:
        self._pending_approval_description = description
        self._pending_approval_args = dict(arguments)
        self._pending_approval_index = 1
        self._pending_approval_total = 1

    def _publish_tool_call_cell(
        self,
        tool_name: str,
        status: str,
        description: str,
        arguments: dict,
        note: str,
    ) -> None:
        call_id = self._direct_tool_call_ids.get(tool_name)
        existing = self._tool_cells.get(call_id or "")
        if call_id is None or (
            status in {"pending", "approved"}
            and existing is not None
            and existing.status
            in {"ok", "error", "denied", "skipped", "interrupted"}
        ):
            call_id = f"direct:{tool_name}:{len(self._tool_order) + 1}"
            self._direct_tool_call_ids[tool_name] = call_id
        self._upsert_tool_cell(
            provider_call_id=call_id,
            tool=tool_name,
            status=status,
            description=description,
            arguments=arguments,
            note=note,
            queue_index=1,
            queue_total=1,
            session_rule_active=tool_name in self._session_allow_tools,
        )

    def _sync_session_allow_tools(self, allowed: set[str]) -> None:
        self._session_allow_tools = set(allowed)

    def _handle_slash_tool_failure(
        self,
        tool_name: str,
        message: str,
        result: dict[str, object] | None = None,
    ) -> None:
        self.query_one(FooterWidget).set_phase(Phase.ERROR)
        self.query_one(FooterWidget).set_hint("Slash command failed")
        if tool_name == "execute_chemsmart_command" and result is not None:
            self._ready_command = None
            self._publish_execute_tool_result(result)
        self.post_error(tool_name, message)

    def _handle_slash_tool_success(
        self,
        tool_name: str,
        arguments: dict[str, object],
        result: object,
    ) -> None:
        footer = self.query_one(FooterWidget)
        footer.set_phase(Phase.FINISHED)
        self._apply_workspace_project_tool_result(tool_name, result)
        if tool_name == "execute_chemsmart_command" and isinstance(
            result, dict
        ):
            self._ready_command = None
            self._publish_execute_tool_result(result)
            footer.set_hint("Command execution completed")
            return
        if tool_name == "wizard_probe" and isinstance(result, dict):
            self._latest_wizard_probe = result
            server_name = str(
                result.get("server_name")
                or arguments.get("server_name")
                or "wizard"
            )
            yaml_text = str(result.get("yaml_text") or "")
            self.post_agent_message(
                f"```yaml\n{yaml_text.rstrip()}\n```",
                title=f"Wizard: {server_name}",
            )
            validation = result.get("validation")
            if isinstance(validation, dict) and not validation.get(
                "ok", False
            ):
                errors = validation.get("errors") or []
                self.post_error(
                    "Wizard validation failed",
                    "\n".join(str(error) for error in errors)
                    or "wizard output did not validate",
                )
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Wizard validation failed")
                return
            footer.set_hint("Wizard YAML ready")
            return
        if tool_name == "wizard_write" and isinstance(result, dict):
            self.post_agent_message(
                f"Wrote `{result.get('written_path')}`.",
                title="Wizard",
            )
            footer.set_hint("Wizard YAML written")
            return
        if tool_name == "write_project_yaml" and isinstance(result, dict):
            project = str(result.get("project_name") or "project")
            program = str(result.get("program") or "gaussian")
            written_path = str(result.get("written_path") or "")
            command_example = (
                f"chemsmart run {program} -p {project} "
                "-f examples/h2o.xyz -c 0 -m 1 opt"
            )
            self.post_agent_message(
                (
                    f"Wrote `{written_path}`.\n\n"
                    "Deterministic use:\n\n"
                    "```bash\n"
                    f"{command_example}\n"
                    "```\n\n"
                    f"This workspace now has `{program}:{project}` loaded. "
                    "The command-synthesis harness will attach and validate "
                    f"this project automatically, and explicit CLI commands "
                    f"can use `-p {project}` from this same workspace."
                ),
                title="Project YAML",
            )
            footer.set_hint(f"Project YAML written: {project}")
            return
        if tool_name == "wizard_refresh" and isinstance(result, dict):
            self.post_agent_message(
                self._wizard_refresh_result_table(result),
                title=f"Wizard refresh: {result.get('server_name') or 'server'}",
            )
            status = str(result.get("status") or "")
            if status == "error":
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Wizard refresh failed")
                return
            if status == "stale":
                footer.set_hint("Wizard refresh preserved stale cache")
                return
            footer.set_hint("Wizard cache ready")
            return
        if tool_name == "wizard_verify" and isinstance(result, dict):
            self.post_agent_message(
                self._wizard_verify_result_table(result),
                title=f"Wizard verify: {result.get('server_name') or 'server'}",
            )
            errors = result.get("errors") or []
            warnings = result.get("warnings") or []
            if errors:
                self.post_error(
                    "Wizard verify failed",
                    "\n".join(str(error) for error in errors),
                )
                footer.set_phase(Phase.ERROR)
                footer.set_hint("Wizard verify found errors")
                return
            if warnings:
                footer.set_hint("Wizard verify completed with warnings")
                return
            footer.set_hint("Wizard verify completed")
            return
        footer.set_hint("Slash command complete")

    def _tool_success_note(self, tool_name: str, result: object) -> str:
        if tool_name == "wizard_probe" and isinstance(result, dict):
            validation = result.get("validation")
            if isinstance(validation, dict) and validation.get("ok", False):
                return "Rendered validated wizard YAML."
            return "Rendered wizard YAML with validation issues."
        if tool_name == "wizard_verify" and isinstance(result, dict):
            if result.get("errors"):
                return "Verified transport wiring and found errors."
            if result.get("warnings"):
                return "Verified transport wiring with warnings."
            return "Verified transport wiring."
        if tool_name == "wizard_refresh" and isinstance(result, dict):
            status = str(result.get("status") or "")
            if status == "error":
                return (
                    "Wizard refresh failed and recorded an error cache entry."
                )
            if status == "stale":
                return "Wizard refresh failed; preserved the last-good stale cache."
            return "Wizard refresh produced a fresh cache entry."
        if tool_name == "wizard_write" and isinstance(result, dict):
            return f"Wrote {result.get('written_path')}"
        if tool_name == "write_project_yaml" and isinstance(result, dict):
            return f"Wrote {result.get('written_path')}"
        return "Completed successfully."

    def _wizard_verify_result_table(self, result: dict[str, object]) -> Table:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Field", style="cyan", no_wrap=True)
        table.add_column("Value")
        table.add_row("server_name", str(result.get("server_name") or "-"))
        table.add_row("host", str(result.get("host") or "-"))
        table.add_row("mode", str(result.get("mode") or "-"))
        table.add_row(
            "would_submit_via",
            str(result.get("would_submit_via") or "-"),
        )
        invocation = result.get("transport_invocation")
        table.add_row(
            "transport_invocation",
            json.dumps(invocation or []),
        )
        warnings = result.get("warnings") or []
        table.add_row(
            "warnings",
            "\n".join(str(item) for item in warnings) if warnings else "-",
        )
        errors = result.get("errors") or []
        table.add_row(
            "errors",
            "\n".join(str(item) for item in errors) if errors else "-",
        )
        return table

    def _wizard_refresh_result_table(self, result: dict[str, object]) -> Table:
        table = Table(show_header=True, box=None, padding=(0, 1))
        table.add_column("Field", style="cyan", no_wrap=True)
        table.add_column("Value")
        table.add_row("cache_path", str(result.get("cache_path") or "-"))
        table.add_row("status", str(result.get("status") or "-"))
        table.add_row("host", str(result.get("host") or "-"))
        table.add_row("mode", str(result.get("mode") or "-"))
        table.add_row("scheduler", str(result.get("scheduler") or "-"))
        table.add_row("probed_at", str(result.get("probed_at") or "-"))
        node_summary = result.get("node_summary")
        if not isinstance(node_summary, dict):
            node_summary = {}
        table.add_row(
            "selected_queue",
            str(node_summary.get("selected_queue") or "-"),
        )
        table.add_row(
            "resources",
            f"cpu={node_summary.get('cpu')} mem_gb={node_summary.get('mem_gb')} gpu={node_summary.get('gpu')}",
        )
        table.add_row(
            "project",
            str(node_summary.get("project") or "-"),
        )
        table.add_row(
            "scratch",
            f"{node_summary.get('scratch_dir') or '-'} (writable={node_summary.get('scratch_writable')})",
        )
        table.add_row("last_error", str(result.get("last_error") or "-"))
        return table

