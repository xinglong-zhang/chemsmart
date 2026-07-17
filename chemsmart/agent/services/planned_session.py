"""Execute a pre-built Plan behind the AgentSession compatibility facade."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from chemsmart.agent.models import CriticVerdict, Plan, SessionState, Step
from chemsmart.agent.services.plan_support import (
    is_project_yaml_workflow,
    render_plan,
)

_RISKY_TOOLS = {"run_local", "submit_hpc"}


class PlannedSessionHost(Protocol):
    state: SessionState | None
    session_dir: Path | None
    decision_log: Any | None
    transport: Any | None

    def _collect_prior_results(
        self, results: list[Any], tool_name: str
    ) -> list[Any]: ...
    def _find_prior_result(
        self, results: list[Any], tool_name: str
    ) -> Any: ...
    def _get_logged_verdict(self) -> CriticVerdict | None: ...
    def _get_logged_summary(self) -> dict[str, Any] | None: ...
    def _execute_step(
        self,
        step_index: int,
        step: Step,
        results: list[Any],
        extra_kwargs: dict[str, Any] | None = None,
    ) -> Any: ...
    def _save_state(self) -> None: ...
    def _preview_submit_step(self, *args: Any, **kwargs: Any) -> Any: ...
    def _critic_call(self, **kwargs: Any) -> CriticVerdict: ...
    def _apply_deterministic_gates(self, **kwargs: Any) -> CriticVerdict: ...
    def _block_reason(self, **kwargs: Any) -> str | None: ...
    def _record_skipped_step(self, *args: Any, **kwargs: Any) -> None: ...
    def _finalize_session(self, **kwargs: Any) -> None: ...


@dataclass
class PlannedRunContext:
    plan: Plan
    results: list[Any]
    dry_run_results: list[dict[str, Any]]
    runtime_result: Any | None
    verdict: CriticVerdict | None
    preview_submit: Any | None = None
    blocked: bool = False
    block_reason: str | None = None
    run_error: Exception | None = None
    paused_for_approval: bool = False
    risky_start: int = field(init=False)

    def __post_init__(self) -> None:
        self.risky_start = len(self.plan.steps)


class PlannedSessionRunner:
    """Run deterministic planner steps and preserve legacy result contracts."""

    def __init__(self, session: PlannedSessionHost) -> None:
        self.session = session

    def run(
        self,
        completed_results: list[Any],
        *,
        dry_submit: bool,
        pause_before_risky: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
        rerender_plan: bool,
    ) -> dict[str, Any]:
        context = self._context(completed_results)
        if not context.plan.steps:
            return self._advisory_result(context, rerender_plan)
        if self._already_completed(context):
            return self._result(context, rerender_plan, include_results=True)

        try:
            self._execute_safe_steps(context)
            self._prepare_submit_preview(context, dry_submit)
            self._ensure_verdict(context, dry_submit)
            self._apply_block_policy(
                context,
                dry_submit=dry_submit,
                allow_remote_unknown=allow_remote_unknown,
                allow_critic_override=allow_critic_override,
            )
            if context.blocked:
                return self._result(context, rerender_plan)
            if self._should_pause(context, pause_before_risky, dry_submit):
                context.paused_for_approval = True
                return self._approval_result(context, rerender_plan)
            self._execute_risky_steps(context, dry_submit)
            return self._result(context, rerender_plan, include_results=True)
        except Exception as exc:
            context.run_error = exc
            raise
        finally:
            if not context.paused_for_approval:
                self.session._finalize_session(
                    verdict=context.verdict,
                    blocked=context.blocked,
                    block_reason=context.block_reason,
                    dry_run_results=context.dry_run_results,
                    run_error=context.run_error,
                )

    def _context(self, completed_results: list[Any]) -> PlannedRunContext:
        session = self.session
        assert session.state is not None
        assert session.state.plan is not None
        assert session.session_dir is not None
        assert session.decision_log is not None
        results = list(completed_results)
        return PlannedRunContext(
            plan=session.state.plan,
            results=results,
            dry_run_results=session._collect_prior_results(
                results, "dry_run_input"
            ),
            runtime_result=session._find_prior_result(
                results, "validate_runtime"
            ),
            verdict=session._get_logged_verdict(),
        )

    def _execute_safe_steps(self, context: PlannedRunContext) -> None:
        session = self.session
        assert session.state is not None
        for step_index in range(
            session.state.current_step_index, len(context.plan.steps)
        ):
            step = context.plan.steps[step_index]
            if step.tool in _RISKY_TOOLS:
                context.risky_start = step_index
                return
            result = session._execute_step(step_index, step, context.results)
            context.results.append(result)
            session.state.current_step_index = step_index + 1
            session._save_state()
            if step.tool == "dry_run_input":
                context.dry_run_results.append(result)
            elif step.tool == "validate_runtime":
                context.runtime_result = result

    def _prepare_submit_preview(
        self,
        context: PlannedRunContext,
        dry_submit: bool,
    ) -> None:
        if context.risky_start >= len(context.plan.steps) or dry_submit:
            return
        step = context.plan.steps[context.risky_start]
        if step.tool == "submit_hpc":
            context.preview_submit = self.session._preview_submit_step(
                context.risky_start,
                step,
                context.results,
                dry_submit=dry_submit,
            )

    def _ensure_verdict(
        self,
        context: PlannedRunContext,
        dry_submit: bool,
    ) -> None:
        if context.verdict is not None or is_project_yaml_workflow(
            context.plan
        ):
            return
        verdict = self.session._critic_call(
            plan=context.plan,
            dry_run_results=context.dry_run_results,
            dry_submit=dry_submit,
        )
        context.verdict = self.session._apply_deterministic_gates(
            plan=context.plan,
            verdict=verdict,
            runtime_result=context.runtime_result,
            dry_run_results=context.dry_run_results,
            preview_submit=context.preview_submit,
            dry_submit=dry_submit,
        )
        assert self.session.decision_log is not None
        self.session.decision_log.write(
            "critic_verdict",
            context.verdict.model_dump(),
            rationale=context.verdict.rationale,
        )

    def _apply_block_policy(
        self,
        context: PlannedRunContext,
        *,
        dry_submit: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
    ) -> None:
        if context.verdict is None:
            return
        context.block_reason = self.session._block_reason(
            verdict=context.verdict,
            dry_submit=dry_submit,
            allow_remote_unknown=allow_remote_unknown,
            allow_critic_override=allow_critic_override,
        )
        context.blocked = context.block_reason is not None

    @staticmethod
    def _should_pause(
        context: PlannedRunContext,
        pause_before_risky: bool,
        dry_submit: bool,
    ) -> bool:
        if not pause_before_risky or context.risky_start >= len(
            context.plan.steps
        ):
            return False
        next_tool = context.plan.steps[context.risky_start].tool
        return not (dry_submit and next_tool == "submit_hpc")

    def _execute_risky_steps(
        self,
        context: PlannedRunContext,
        dry_submit: bool,
    ) -> None:
        session = self.session
        assert session.state is not None
        for step_index in range(context.risky_start, len(context.plan.steps)):
            if step_index < session.state.current_step_index:
                continue
            step = context.plan.steps[step_index]
            if step.tool == "submit_hpc" and dry_submit:
                session._record_skipped_step(
                    step_index,
                    step,
                    reason="dry-submit skips remote submission",
                )
                context.preview_submit = {
                    "skipped": True,
                    "skip_reason": "dry-submit skips remote submission",
                }
                return
            extra_kwargs = self._execution_kwargs(step, dry_submit)
            result = session._execute_step(
                step_index,
                step,
                context.results,
                extra_kwargs=extra_kwargs,
            )
            context.results.append(result)
            session.state.current_step_index = step_index + 1
            session._save_state()

    def _execution_kwargs(
        self,
        step: Step,
        dry_submit: bool,
    ) -> dict[str, Any]:
        if step.tool != "submit_hpc":
            return {}
        kwargs: dict[str, Any] = {"execute": not dry_submit}
        if self.session.transport is not None:
            kwargs["transport"] = self.session.transport
        return kwargs

    def _already_completed(self, context: PlannedRunContext) -> bool:
        assert self.session.state is not None
        return (
            self.session.state.current_step_index >= len(context.plan.steps)
            and self.session._get_logged_summary() is not None
        )

    def _advisory_result(
        self,
        context: PlannedRunContext,
        rerender_plan: bool,
    ) -> dict[str, Any]:
        is_chitchat = context.plan.is_chitchat()
        self.session._finalize_session(
            verdict=None,
            blocked=False,
            block_reason=None,
            dry_run_results=context.dry_run_results,
            advisory_only=True,
            is_chitchat=is_chitchat,
            rationale=context.plan.rationale,
        )
        result = self._result(context, rerender_plan, include_results=True)
        result.update(
            {
                "dry_run_result": None,
                "advisory_only": True,
                "is_chitchat": is_chitchat,
            }
        )
        return result

    def _approval_result(
        self,
        context: PlannedRunContext,
        rerender_plan: bool,
    ) -> dict[str, Any]:
        result = self._result(context, rerender_plan)
        result["pending_approval"] = True
        result["next_risky_tool"] = context.plan.steps[
            context.risky_start
        ].tool
        return result

    def _result(
        self,
        context: PlannedRunContext,
        rerender_plan: bool,
        *,
        include_results: bool = False,
    ) -> dict[str, Any]:
        assert self.session.state is not None
        assert self.session.session_dir is not None
        result = {
            "session_id": self.session.state.session_id,
            "session_dir": str(self.session.session_dir),
            "plan": context.plan,
            "plan_text": render_plan(context.plan) if rerender_plan else None,
            "critic_verdict": context.verdict,
            "completed_steps": self.session.state.current_step_index,
            "blocked": context.blocked,
            "dry_run_result": _primary_dry_run_result(context.dry_run_results),
            "dry_run_results": context.dry_run_results,
            "runtime_result": context.runtime_result,
            "preview_submit": context.preview_submit,
        }
        if include_results:
            result["results"] = context.results
        return result


def _primary_dry_run_result(
    dry_run_results: list[dict[str, Any]],
) -> dict[str, Any] | None:
    return dry_run_results[0] if dry_run_results else None


__all__ = ["PlannedSessionRunner"]
