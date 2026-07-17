"""Background workers that call agent, synthesis, and project services."""

from __future__ import annotations

from textual import work

from chemsmart.agent.core import AgentSession
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.synthesis import SynthesisSession
from chemsmart.agent.tui.chat_helpers import (
    _ExecuteOverrideRegistry,
    _extract_execute_context,
    _format_synthesis_exception,
    _new_synthesis_artifact_dir,
    _provider_model_label,
    _provider_type_label,
    _write_synthesis_artifact,
)


class SessionWorkersMixin:
    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-unified",
    )
    def run_unified_session(self, request: str) -> dict[str, object]:
        self.active_resume_id = None
        if (
            self.active_agent_session is None
            or self.active_agent_session.stage_prompt != "unified_agent.md"
        ):
            self.active_agent_session = AgentSession(
                session_root=str(self.session_root),
                stage_prompt="unified_agent.md",
                runtime_v2=self.runtime_v2,
            )
        policy = self._permission_policy(prompt_risky=True)
        result = self.active_agent_session.run_loop(
            request,
            policy=policy,
            approver=lambda req: self._await_approval(req),
        )
        result["session_allow_tools"] = sorted(policy.session_allow)
        return result

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-ask",
    )
    def run_synthesis_session(self, request: str) -> dict[str, object]:
        provider_type = _provider_type_label(self._active_provider_config)
        provider_model = _provider_model_label(self._active_provider_config)
        artifact_dir = _new_synthesis_artifact_dir(
            self.session_root,
            provider_type=provider_type,
        )
        try:
            if self.active_synthesis_session is None:
                self.active_synthesis_session = SynthesisSession()
            result = self.active_synthesis_session.prepare_command(request)
            semantic = self.active_synthesis_session._last_semantic_result
            payload = {
                "request": request,
                "provider_type": provider_type,
                "provider_model": provider_model,
                "synthesis": result,
                "raw_response": self.active_synthesis_session._last_raw_response,
                "semantic_result": (
                    semantic.to_dict() if semantic is not None else None
                ),
            }
            _write_synthesis_artifact(artifact_dir, payload)
        except Exception as exc:
            payload = {
                "synthesis": {
                    "status": "infeasible",
                    "command": "",
                    "explanation": _format_synthesis_exception(exc),
                    "confidence": "low",
                    "missing_info": [],
                    "alternatives": [],
                },
                "raw_response": "",
                "semantic_result": None,
            }
            _write_synthesis_artifact(
                artifact_dir,
                {
                    "request": request,
                    "provider_type": provider_type,
                    "provider_model": provider_model,
                    **payload,
                },
            )
            return {
                **payload,
                "provider_type": provider_type,
                "provider_model": provider_model,
                "synthesis_artifact_dir": str(artifact_dir),
            }
        return {
            "synthesis": result,
            "raw_response": self.active_synthesis_session._last_raw_response,
            "semantic_result": (
                semantic.to_dict() if semantic is not None else None
            ),
            "provider_type": provider_type,
            "provider_model": provider_model,
            "synthesis_artifact_dir": str(artifact_dir),
        }

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="slash-tool",
        name="project-yaml-write",
    )
    def run_project_yaml_write(
        self,
        arguments: dict[str, object],
    ) -> dict[str, object]:
        tool_name = "write_project_yaml"
        registry = ToolRegistry.default()
        normalized_args = registry.normalize_args(tool_name, arguments)
        description = registry.describe_tool(tool_name)
        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "approved",
            description,
            normalized_args,
            "Confirmed in the workspace project write dialog.",
        )
        result = registry.call(tool_name, normalized_args)
        if isinstance(result, dict) and result.get("ok") is False:
            message = str(result.get("error") or "Project YAML write failed.")
            self.app.call_from_thread(
                self._publish_tool_call_cell,
                tool_name,
                "error",
                description,
                normalized_args,
                message,
            )
            self.app.call_from_thread(
                self._handle_slash_tool_failure,
                tool_name,
                message,
            )
            return {"tool": tool_name, "status": "error", "result": result}

        self.app.call_from_thread(
            self._publish_tool_call_cell,
            tool_name,
            "ok",
            description,
            normalized_args,
            self._tool_success_note(tool_name, result),
        )
        self.app.call_from_thread(
            self._handle_slash_tool_success,
            tool_name,
            normalized_args,
            result,
        )
        return {"tool": tool_name, "status": "ok", "result": result}

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
    ) -> dict[str, object]:
        self.active_resume_id = session_id
        self.active_agent_session = AgentSession.load(
            session_id,
            session_root=str(self.session_root),
            cwd_override=cwd_override,
            runtime_v2=self.runtime_v2,
        )
        policy = self._permission_policy()
        request = self.active_agent_session.state.request or "Continue."
        result = self.active_agent_session.run_loop(
            request,
            policy=policy,
            approver=lambda req: self._await_approval(req),
        )
        result["session_allow_tools"] = sorted(policy.session_allow)
        return result

    @work(
        thread=True,
        exclusive=True,
        exit_on_error=False,
        group="agent-session",
        name="agent-execute",
    )
    def execute_agent_session(self, session_id: str) -> dict[str, object]:
        from chemsmart.agent.transport import SshQsubTransport

        self.active_resume_id = session_id
        session_dir = self.session_root / session_id
        if self.active_agent_session is None or str(
            self.active_agent_session.session_dir
        ) != str(session_dir):
            self.active_agent_session = AgentSession.load(
                session_id,
                session_root=str(self.session_root),
                runtime_v2=self.runtime_v2,
            )
        session = self.active_agent_session

        ctx = _extract_execute_context(session_dir)
        job_handle_id = ctx.get("job_handle_id")
        if job_handle_id:
            parts = [
                "The dry-run has been reviewed and approved by the user.",
                f"The job is stored as handle '{job_handle_id}'.",
            ]
            if ctx.get("inputfile"):
                parts.append(f"The input file is at '{ctx['inputfile']}'.")
            if ctx.get("server_name"):
                parts.append(f"Server name: '{ctx['server_name']}'.")
            parts.append(f"Please call submit_hpc with job='{job_handle_id}'.")
            execute_request = " ".join(parts)
        else:
            execute_request = (
                "The dry-run and all pre-flight checks have been reviewed "
                "and approved by the user. Please submit the job to the "
                "HPC scheduler now."
            )

        original_registry = session.registry
        session.registry = _ExecuteOverrideRegistry(
            original_registry, SshQsubTransport()
        )
        policy = PermissionPolicy(
            mode=PermissionMode.DRIVING,
            yolo=False,
            session_allow=set(),
            driving_denylist=set(),
        )
        try:
            result = session.run_loop(
                execute_request,
                policy=policy,
                approver=lambda req: ApprovalDecision.ALLOW_ONCE,
            )
        finally:
            session.registry = original_registry
        result["session_allow_tools"] = []
        return result

    def _load_agent_session(
        self,
        session_id: str,
        *,
        cwd_override: str | None = None,
    ) -> AgentSession:
        self.active_resume_id = session_id
        self.active_agent_session = AgentSession.load(
            session_id,
            session_root=str(self.session_root),
            cwd_override=cwd_override,
            runtime_v2=self.runtime_v2,
        )
        return self.active_agent_session

