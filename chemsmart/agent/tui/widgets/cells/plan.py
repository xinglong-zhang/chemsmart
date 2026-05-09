"""Plan/workflow transcript cell."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rich.console import Group
from rich.text import Text

from chemsmart.agent.core import Plan, Step, _restore_json_result
from chemsmart.io.molecules.structure import Molecule

from .base import BaseCell, no_data_text


@dataclass(slots=True)
class _WorkflowRow:
    index: int
    tool: str
    args: dict[str, Any]
    detail: str
    status: str = "pending"


_STATUS_ICON = {
    "done": "✓",
    "running": "⠋",
    "failed": "✖",
    "pending": "•",
}


class PlanCell(BaseCell):
    def __init__(self, content: Plan | str) -> None:
        self.plan = content if isinstance(content, Plan) else None
        self.plan_text = content if isinstance(content, str) else None
        self._rows = (
            [
                _row_from_step(index, step)
                for index, step in enumerate(
                    self.plan.steps,
                    start=1,
                )
            ]
            if self.plan is not None
            else []
        )
        super().__init__(
            self._build_renderable(),
            title="Workflow" if self.plan is not None else "Plan",
            classes="plan-cell",
        )

    def mark_started(
        self,
        step_index: int,
        tool: str,
        args: dict[str, Any] | None = None,
    ) -> None:
        row = self._find_row(step_index, tool)
        if row is None:
            return
        row.status = "running"
        row.detail = _running_detail(tool, args or row.args)
        self.update(self._build_renderable())

    def mark_completed(
        self,
        step_index: int,
        tool: str,
        payload: Any | None = None,
    ) -> None:
        row = self._find_row(step_index, tool)
        if row is None:
            return
        row.status = "done"
        row.detail = _completed_detail(tool, row.args, payload)
        self.update(self._build_renderable())

    def mark_failed(
        self,
        step_index: int,
        tool: str,
        message: str,
    ) -> None:
        row = self._find_row(step_index, tool)
        if row is None:
            return
        row.status = "failed"
        row.detail = message.strip() or "실행이 중단되었습니다."
        self.update(self._build_renderable())

    def _find_row(self, step_index: int, tool: str) -> _WorkflowRow | None:
        for row in self._rows:
            if row.index == step_index and row.tool == tool:
                return row
        return None

    def _build_renderable(self):
        if self.plan is None:
            text = self.plan_text or ""
            return Text(text) if text.strip() else no_data_text()
        lines = []
        for row in self._rows:
            icon = _STATUS_ICON.get(row.status, "•")
            line = Text()
            style = _status_style(row.status)
            line.append(f"{icon} ", style=style)
            line.append(f"{row.index}. {row.tool}", style="bold")
            if row.detail:
                line.append("  ", style="dim")
                line.append(row.detail, style=style)
            lines.append(line)
        return Group(*lines)


def _row_from_step(index: int, step: Step) -> _WorkflowRow:
    return _WorkflowRow(
        index=index,
        tool=step.tool,
        args=dict(step.args or {}),
        detail=_pending_detail(step.tool, dict(step.args or {})),
    )


def _status_style(status: str) -> str:
    if status == "done":
        return "success"
    if status == "running":
        return "accent"
    if status == "failed":
        return "error"
    return "dim"


def _pending_detail(tool: str, args: dict[str, Any]) -> str:
    if tool == "build_molecule":
        path = str(args.get("filepath") or "")
        return f"{path} 대기 중" if path else "구조 파일 대기 중"
    if tool == "dry_run_input":
        return "입력 파일 렌더링 대기 중"
    if tool == "validate_runtime":
        return "런타임 점검 대기 중"
    if tool == "run_local":
        return "실행 대기 중"
    if tool == "submit_hpc":
        return "제출 미리보기 대기 중"
    return "대기 중"


def _running_detail(tool: str, args: dict[str, Any]) -> str:
    if tool == "build_molecule":
        path = str(args.get("filepath") or "")
        return f"{path} 확인 중…" if path else "구조 파일 확인 중…"
    if tool == "build_gaussian_settings":
        return "Gaussian 설정 정리 중…"
    if tool == "build_job":
        return "계산 작업 조립 중…"
    if tool == "dry_run_input":
        return "입력 파일 렌더링 중…"
    if tool == "validate_runtime":
        return "런타임 점검 중…"
    if tool == "run_local":
        return "로컬 실행 시작 중…"
    if tool == "submit_hpc":
        return "제출 미리보기 생성 중…"
    return "진행 중…"


def _completed_detail(
    tool: str,
    args: dict[str, Any],
    payload: Any | None,
) -> str:
    restored = _restore_json_result(payload)
    if tool == "build_molecule" and isinstance(restored, Molecule):
        formula = (
            restored.formula() if hasattr(restored, "formula") else "Molecule"
        )
        count = len(restored.symbols or [])
        charge = restored.charge if restored.charge is not None else 0
        multiplicity = (
            restored.multiplicity if restored.multiplicity is not None else 1
        )
        source = Path(str(args.get("filepath") or "")).as_posix()
        details = f"{formula} · {count} atoms · q={charge} mult={multiplicity}"
        return f"{details} · {source}" if source else details
    if tool == "build_gaussian_settings" and isinstance(restored, dict):
        functional = str(restored.get("functional") or "").upper()
        basis = str(restored.get("basis") or "").upper()
        return f"{functional}/{basis}".strip("/")
    if tool == "build_gaussian_settings":
        functional = str(getattr(restored, "functional", "") or "").upper()
        basis = str(getattr(restored, "basis", "") or "").upper()
        return f"{functional}/{basis}".strip("/")
    if tool == "build_job":
        job = restored
        kind = str(args.get("kind") or "")
        label = str(getattr(job, "label", "") or args.get("label") or "")
        return " · ".join(part for part in (kind, label) if part)
    if tool == "dry_run_input" and isinstance(restored, dict):
        inputfile = str(restored.get("inputfile") or "")
        return Path(inputfile).name if inputfile else "입력 파일 준비됨"
    if tool == "validate_runtime" and isinstance(restored, dict):
        return _runtime_summary(restored)
    if tool == "run_local" and isinstance(restored, dict):
        if restored.get("ok"):
            return "로컬 실행 완료"
        return f"returncode {restored.get('returncode', '?')}"
    if tool == "submit_hpc" and isinstance(restored, dict):
        if restored.get("job_id"):
            return f"job {restored['job_id']}"
        return "제출 미리보기 준비됨"
    return "완료"


def _runtime_summary(validation: dict[str, Any]) -> str:
    remote_unknown = validation.get("remote_unknown") or []
    local_issues = validation.get("local_issues") or []
    if local_issues:
        return f"로컬 이슈 {len(local_issues)}개"
    if remote_unknown:
        return f"로컬 OK / 원격 정보 {len(remote_unknown)}개 필요"
    return "로컬/원격 실행 준비됨"
