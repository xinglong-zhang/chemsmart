"""Validated synthesis state and transcript presentation."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.tui.chat_helpers import (
    _command_details_text,
    _decision_trace_dict,
    _final_answer_text,
    _final_command_text,
    _format_semantic_result,
    _ready_command_hint,
)
from chemsmart.agent.tui.chat_models import ReadyCommand
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CommandInterpretationCell,
    DecisionTraceCell,
    FinalAnswerCell,
    SynthesisTraceCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    FilePickerOverlay,
    TextPromptOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript
from chemsmart.settings.workspace_project import workspace_project_path


class SynthesisPresentationMixin:
    """Promote gated commands and present synthesis outcomes."""

    def _can_approve_or_execute(self, action: str) -> bool:
        pending = self._pending_tool_request
        if self._pending_approval and pending is not None:
            if pending.name == "execute_chemsmart_command":
                pending_command = str(pending.arguments.get("command") or "")
                return parse_model_command(pending_command).action == action
            return pending.name == (
                "run_local" if action == "run" else "submit_hpc"
            )
        return bool(
            self._ready_command is not None
            and self._ready_command.action == action
        )

    def _remember_ready_command(
        self,
        *,
        command: str,
        semantic: dict[str, object] | None,
        intent: dict[str, object] | None,
        source: str,
    ) -> bool:
        self._ready_command = None
        parsed = parse_model_command(command)
        if parsed.parse_error or parsed.action not in {"run", "sub"}:
            return False
        if not isinstance(semantic, dict):
            return False
        semantic_verdict = str(semantic.get("verdict") or "").lower()
        if semantic_verdict not in {"ok", "warn"}:
            return False
        generated_inputs = semantic.get("generated_inputs")
        if not isinstance(generated_inputs, list) or not generated_inputs:
            return False
        if not isinstance(intent, dict):
            return False
        intent_verdict = str(intent.get("verdict") or "").lower()
        if intent_verdict not in {"ok", "warn"}:
            return False

        project_path: Path | None = None
        project_hash: str | None = None
        if parsed.program in {"gaussian", "orca"}:
            if not parsed.project:
                return False
            project_path = workspace_project_path(
                parsed.project,
                parsed.program,
                cwd=Path.cwd(),
            )
            if not project_path.is_file():
                return False
            project_hash = sha256(project_path.read_bytes()).hexdigest()

        self._ready_command = ReadyCommand(
            command=command,
            action=parsed.action,
            workspace=Path.cwd().resolve(),
            semantic_verdict=semantic_verdict,
            intent_verdict=intent_verdict,
            project_path=project_path,
            project_sha256=project_hash,
            source=source,
        )
        return True

    def _ready_command_problem(self, expected_action: str) -> str | None:
        ready = self._ready_command
        if ready is None:
            return "No semantic- and intent-validated command is ready."
        if ready.action != expected_action:
            expected = (
                "chemsmart run"
                if expected_action == "run"
                else "chemsmart sub"
            )
            return f"The validated command is not a `{expected}` command."
        if ready.workspace != Path.cwd().resolve():
            return (
                "The workspace changed after validation. Regenerate the command "
                "from the current workspace before execution."
            )
        if ready.project_path is not None:
            if not ready.project_path.is_file():
                return "The validated workspace project YAML no longer exists."
            current_hash = sha256(ready.project_path.read_bytes()).hexdigest()
            if current_hash != ready.project_sha256:
                return (
                    "The workspace project YAML changed after validation. "
                    "Regenerate the command so the generated-input evidence "
                    "matches the current project."
                )
        return None

    def _publish_synthesis_result(self, result: dict[str, object]) -> None:
        synthesis = result.get("synthesis")
        if not isinstance(synthesis, dict):
            self.post_error("Synthesis failed", "Provider returned no result.")
            self.query_one(FooterWidget).set_phase(Phase.ERROR)
            self.query_one(FooterWidget).set_hint("Synthesis failed")
            return

        status = str(synthesis.get("status") or "")
        self._waiting_for_user = status == "needs_clarification"
        footer = self.query_one(FooterWidget)
        provider_type = str(result.get("provider_type") or "offline")
        provider_model = str(result.get("provider_model") or "auto")
        artifact_dir = str(result.get("synthesis_artifact_dir") or "")
        semantic = result.get("semantic_result")
        semantic_dict = semantic if isinstance(semantic, dict) else None
        if status == "ready":
            command = str(synthesis.get("command") or "")
            explanation = str(synthesis.get("explanation") or "")
            confidence = str(synthesis.get("confidence") or "low")
            project = str(synthesis.get("project") or "")
            intent = synthesis.get("intent_assertion") or synthesis.get(
                "intent"
            )
            intent_dict = intent if isinstance(intent, dict) else None
            command_is_executable = self._remember_ready_command(
                command=command,
                semantic=semantic_dict,
                intent=intent_dict,
                source="local_synthesis",
            )
            transcript = self.query_one(Transcript)
            transcript.add_cell(
                SynthesisTraceCell(
                    provider_type=provider_type,
                    model=provider_model,
                    mode=self._interaction_mode,
                    status=status,
                    command=command,
                    semantic=semantic_dict,
                    decision_trace=_decision_trace_dict(synthesis),
                    artifact_dir=artifact_dir,
                )
            )
            self._publish_decision_trace(synthesis)
            transcript.add_cell(
                CommandInterpretationCell(
                    parse_model_command(command),
                    expanded=True,
                )
            )
            transcript.add_cell(
                AgentMessageCell(
                    _command_details_text(
                        explanation=(
                            explanation or "Prepared chemsmart command."
                        ),
                        confidence=confidence,
                        project=project,
                    ),
                    title="Command details",
                )
            )
            transcript.add_cell(
                FinalAnswerCell(
                    _final_command_text(command=command),
                    title="Final Command",
                ),
            )
            footer.set_phase(Phase.FINISHED)
            footer.set_hint(
                _ready_command_hint(self._ready_command)
                if command_is_executable
                else "Command shown, but execution evidence is incomplete"
            )
            transcript.collapse_tool_chain(self._active_turn_id)
            return

        if status == "informational":
            command = str(synthesis.get("command") or "")
            explanation = str(synthesis.get("explanation") or "")
            action = str(synthesis.get("action") or "explain_command")
            title = {
                "explain_command": "Command Explanation",
                "critique_command": "Command Critic",
                "repair_command": "Command Repair",
            }.get(action, "Command Analysis")
            transcript = self.query_one(Transcript)
            transcript.add_cell(
                SynthesisTraceCell(
                    provider_type=provider_type,
                    model=provider_model,
                    mode=self._interaction_mode,
                    status=status,
                    command=command,
                    semantic=semantic_dict,
                    decision_trace=_decision_trace_dict(synthesis),
                    artifact_dir=artifact_dir,
                )
            )
            self._publish_decision_trace(synthesis)
            if command:
                transcript.add_cell(
                    CommandInterpretationCell(
                        parse_model_command(command),
                        expanded=True,
                    )
                )
            transcript.add_cell(
                FinalAnswerCell(
                    _final_answer_text(
                        command=command,
                        explanation=(
                            explanation or "No explanation was generated."
                        ),
                    ),
                    title=title,
                )
            )
            footer.set_phase(Phase.FINISHED)
            footer.set_hint("Command analysis ready")
            transcript.collapse_tool_chain(self._active_turn_id)
            return

        if status == "needs_clarification":
            missing = synthesis.get("missing_info") or []
            if not isinstance(missing, list):
                missing = [str(missing)]
            lines = "\n".join(f"- {item}" for item in missing) or "- details"
            self.query_one(Transcript).add_cell(
                SynthesisTraceCell(
                    provider_type=provider_type,
                    model=provider_model,
                    mode=self._interaction_mode,
                    status=status,
                    command=str(synthesis.get("command") or ""),
                    semantic=semantic_dict,
                    decision_trace=_decision_trace_dict(synthesis),
                    artifact_dir=artifact_dir,
                )
            )
            self._publish_decision_trace(synthesis)
            self.query_one(Transcript).add_cell(
                FinalAnswerCell(
                    (
                        "I need more information before making a command:\n\n"
                        f"{lines}"
                    ),
                    title="Clarification",
                )
            )
            footer.set_phase(Phase.WAITING_USER)
            footer.set_hint("Clarification needed")
            self.query_one(Transcript).collapse_tool_chain(
                self._active_turn_id
            )
            return

        explanation = str(
            synthesis.get("explanation") or "No executable command was made."
        )
        semantic_text = _format_semantic_result(semantic_dict)
        self.query_one(Transcript).add_cell(
            SynthesisTraceCell(
                provider_type=provider_type,
                model=provider_model,
                mode=self._interaction_mode,
                status=status,
                command=str(synthesis.get("command") or ""),
                semantic=semantic_dict,
                decision_trace=_decision_trace_dict(synthesis),
                artifact_dir=artifact_dir,
            )
        )
        self._publish_decision_trace(synthesis)
        self.query_one(Transcript).add_cell(
            FinalAnswerCell(
                f"{explanation}{semantic_text}",
                title="Final Status",
            )
        )
        footer.set_phase(Phase.FINISHED)
        footer.set_hint("No command generated")
        self.query_one(Transcript).collapse_tool_chain(self._active_turn_id)

    def _publish_decision_trace(self, synthesis: dict[str, object]) -> None:
        trace = synthesis.get("decision_trace")
        if isinstance(trace, dict) and trace:
            self.query_one(Transcript).add_cell(DecisionTraceCell(trace))

    def open_file_picker(self) -> None:
        self.app.push_screen(
            FilePickerOverlay(Path.cwd()), self._handle_file_pick
        )

    def edit_method_from_cell(self, recommendation: dict) -> None:
        method = recommendation.get("functional") or "manual"
        basis = recommendation.get("basis") or "manual"
        self.app.push_screen(
            TextPromptOverlay(
                title="Revise method",
                prompt=(
                    "Describe the corrected method/basis to use. "
                    f"Current: {method}/{basis}"
                ),
            ),
            self._handle_method_revision,
        )

    def _handle_file_pick(self, value: str | None) -> None:
        if not value:
            return
        self.query_one(Composer).insert_file_reference(value)

    def _handle_method_revision(self, value: str | None) -> None:
        if not value or not self._current_request:
            return
        self.start_request(self._corrected_request(value))
