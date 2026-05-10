from __future__ import annotations

import asyncio

from textual.widgets import Static

from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.popups import build_approval_overlay


def test_approval_overlay_accepts_arbitrary_tool_and_args(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        results = []
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                build_approval_overlay(
                    tool_name="build_orca_settings",
                    description="Build validated ORCA job settings.",
                    arguments={"functional": "wb97x-d", "charge": 0},
                    session_rule_active=False,
                    queue_index=1,
                    queue_total=2,
                ),
                results.append,
            )
            await pilot.pause()
            summary = app.screen.query_one("#approval-summary", Static)
            text = str(summary.renderable)
            assert "build_orca_settings" in text
            assert '"charge": 0' in text
            assert "1 of 2" in text
            await pilot.press("y")
            await pilot.pause()

        assert results[0].to_decision() == ApprovalDecision.ALLOW_ONCE

    asyncio.run(scenario())


def test_approval_overlay_session_and_deny_shortcuts_map_to_decisions(
    tmp_path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        session_results = []
        deny_results = []
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                build_approval_overlay(
                    tool_name="dry_run_input",
                    description="Render an input file.",
                    arguments={"job": "job_1"},
                    session_rule_active=True,
                ),
                session_results.append,
            )
            await pilot.pause()
            await pilot.press("s")
            await pilot.pause()

            app.push_screen(
                build_approval_overlay(
                    tool_name="run_local",
                    description="Execute a local run.",
                    arguments={"job": "job_1"},
                    session_rule_active=False,
                ),
                deny_results.append,
            )
            await pilot.pause()
            await pilot.press("n")
            await pilot.pause()

        assert (
            session_results[0].to_decision() == ApprovalDecision.ALLOW_SESSION
        )
        assert deny_results[0].to_decision() == ApprovalDecision.DENY

    asyncio.run(scenario())


def test_permission_policy_session_rule_persists_after_allow_session():
    policy = PermissionPolicy(mode=PermissionMode.PERMISSION)

    policy.record("dry_run_input", ApprovalDecision.ALLOW_SESSION)

    assert policy.session_allow == {"dry_run_input"}
