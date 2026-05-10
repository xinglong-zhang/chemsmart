from __future__ import annotations

import asyncio
from threading import Event

from chemsmart.agent.permissions import ApprovalDecision, PermissionMode
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.tui.app import ChemsmartTuiApp


def _pending_request(name: str) -> ToolRequest:
    return ToolRequest(
        request_id=f"openai:{name}",
        provider="openai",
        provider_call_id=f"{name}-1",
        name=name,
        arguments_json="{}",
        arguments={},
        raw={},
    )


def test_permission_slash_commands_route_and_backcompat(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        resolved: list[ApprovalDecision] = []
        requested: list[str] = []

        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._resolve_pending_approval = (
                lambda decision: resolved.append(decision) or True
            )
            screen._request_approval = requested.append

            screen._handle_slash_command("/allow")
            screen._handle_slash_command("/allow-session")
            screen._handle_slash_command("/deny")
            assert resolved == [
                ApprovalDecision.ALLOW_ONCE,
                ApprovalDecision.ALLOW_SESSION,
                ApprovalDecision.DENY,
            ]

            screen._handle_slash_command("/permissions permission")
            assert screen._permission_mode == PermissionMode.PERMISSION
            screen._handle_slash_command("/permissions driving")
            assert screen._permission_mode == PermissionMode.DRIVING

            screen._handle_slash_command("/yolo on")
            assert screen._yolo_enabled is True
            screen._handle_slash_command("/yolo off")
            assert screen._yolo_enabled is False

            screen._handle_slash_command("/run")
            screen._handle_slash_command("/submit")
            assert requested == ["run_local", "submit_hpc"]

            screen._handle_slash_command("/run yes")
            screen._handle_slash_command("/submit session")
            assert resolved[-2:] == [
                ApprovalDecision.ALLOW_ONCE,
                ApprovalDecision.ALLOW_SESSION,
            ]

            screen._pending_approval = True
            screen._approval_waiter = Event()
            screen._pending_tool_request = _pending_request("build_molecule")
            screen._handle_plain_approval_alias("y")
            screen._pending_approval = True
            screen._approval_waiter = Event()
            screen._pending_tool_request = _pending_request("dry_run_input")
            screen._handle_plain_approval_alias("session")
            screen._pending_approval = True
            screen._approval_waiter = Event()
            screen._pending_tool_request = _pending_request("run_local")
            screen._handle_plain_approval_alias("n")
            assert resolved[-3:] == [
                ApprovalDecision.ALLOW_ONCE,
                ApprovalDecision.ALLOW_SESSION,
                ApprovalDecision.DENY,
            ]

    asyncio.run(scenario())


def test_wizard_slash_commands_route_probe_and_write(tmp_path):
    class DummyWorker:
        is_finished = True

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        calls: list[tuple[str, dict]] = []

        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen.run_slash_tool_request = (  # type: ignore[method-assign]
                lambda tool, args: calls.append((tool, args)) or DummyWorker()
            )

            screen._handle_slash_command("/wizard perlmutter login.cluster")
            assert calls == [
                (
                    "wizard_probe",
                    {
                        "server_name": "perlmutter",
                        "ssh_host_hint": "login.cluster",
                    },
                )
            ]

            screen._latest_wizard_probe = {
                "server_name": "perlmutter",
                "yaml_text": "SERVER:\n  SCHEDULER: SLURM\n",
                "validation": {"ok": True, "errors": []},
            }
            screen._handle_slash_command("/wizard-write overwrite")

            assert calls[-1] == (
                "wizard_write",
                {
                    "server_name": "perlmutter",
                    "yaml_text": "SERVER:\n  SCHEDULER: SLURM\n",
                    "overwrite": True,
                },
            )

    asyncio.run(scenario())
