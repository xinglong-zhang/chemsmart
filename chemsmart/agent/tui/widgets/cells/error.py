"""Error transcript cell."""

from __future__ import annotations

import json
from dataclasses import dataclass

from rich.console import Group
from rich.text import Text
from textual.binding import Binding

from .base import BaseCell


@dataclass(slots=True)
class _NormalizedError:
    title: str
    message: str
    notes: list[str]


class ErrorCell(BaseCell):
    BINDINGS = [Binding("e", "toggle_details", "Expand", show=False)]

    def __init__(
        self,
        title: str,
        message: str,
        details: dict | str | None = None,
    ) -> None:
        self.error_title = title or "Error"
        self.message = message or "no data"
        self.details = details
        self.show_details = False
        super().__init__(
            self._build_renderable(),
            title="Error",
            classes="error-cell",
        )

    def action_toggle_details(self) -> None:
        if not self._details_text():
            return
        self.show_details = not self.show_details
        self.update(self._build_renderable())

    def _details_text(self) -> str:
        if self.details is None:
            return ""
        if isinstance(self.details, str):
            return self.details.strip()
        payload = dict(self.details)
        stage = payload.pop("stage", None) or payload.pop("tool", None)
        parts = []
        if stage:
            parts.append(f"Context: {stage}")
        stack = payload.pop("traceback", None) or payload.pop("stack", None)
        if payload:
            parts.append(json.dumps(payload, indent=2, sort_keys=True))
        if stack:
            parts.append(str(stack))
        return "\n\n".join(part for part in parts if part).strip()

    def _build_renderable(self):
        normalized = _normalize_error(
            self.error_title,
            self.message,
            self.details,
        )
        summary = Text.assemble(
            ("✖ ", "bold error"),
            (normalized.title, "bold error"),
        )
        lines = [summary, Text(normalized.message)]
        lines.extend(Text(note, style="dim") for note in normalized.notes)

        details = self._details_text()
        if not details:
            return Group(*lines)

        hint = "  [e expand]" if not self.show_details else "  [e collapse]"
        lines[0].append(hint, style="dim")
        if not self.show_details:
            return Group(*lines)
        lines.extend([Text(""), Text(details, style="dim")])
        return Group(*lines)


def _normalize_error(
    title: str,
    message: str,
    details: dict | str | None,
) -> _NormalizedError:
    payload = details if isinstance(details, dict) else {}
    lowered_message = message.lower()
    if payload.get("tool") == "build_molecule" and (
        "file" in title.lower()
        or "no such file" in lowered_message
        or "not found" in lowered_message
    ):
        return _NormalizedError(
            title="입력 파일을 찾을 수 없습니다",
            message=(
                "요청에 적힌 경로가 현재 작업 폴더 기준으로 존재하지 "
                "않습니다."
            ),
            notes=["예: examples/h2o.xyz 처럼 상대경로를 다시 확인하세요."],
        )
    if title == "Resume failed" and message.startswith("Unknown session id:"):
        return _NormalizedError(
            title="세션을 찾지 못했습니다",
            message=("해당 session id가 현재 저장소에서 보이지 않습니다."),
            notes=[
                (
                    "/sessions 로 최근 세션을 확인하거나, 다른 작업 "
                    "디렉터리에서 생성된 세션인지 확인하세요."
                )
            ],
        )
    if title == "Resume blocked" and isinstance(payload, dict):
        notes = [
            (
                "이 세션은 다른 폴더에서 만들어졌습니다. 경로 해석 "
                "결과가 달라질 수 있습니다."
            )
        ]
        recorded_cwd = payload.get("recorded_cwd")
        current_cwd = payload.get("current_cwd")
        if recorded_cwd and current_cwd:
            notes.extend(
                [
                    f"recorded: {recorded_cwd}",
                    f"current:  {current_cwd}",
                ]
            )
        return _NormalizedError(
            title="작업 디렉터리가 다릅니다",
            message=(
                "이 세션은 다른 폴더에서 만들어졌습니다. 경로 해석 결과가 "
                "달라질 수 있습니다."
            ),
            notes=notes,
        )
    lowered = f"{title} {message}".lower()
    if "registry" in lowered or "schema" in lowered:
        return _NormalizedError(
            title="도구 구성이 예상과 다릅니다",
            message=("등록된 도구 수나 스키마가 이전 세션과 달라졌습니다."),
            notes=[
                (
                    "이전 transcript는 읽을 수 있지만, 재실행 결과는 "
                    "달라질 수 있습니다."
                )
            ],
        )
    return _NormalizedError(
        title=title or "Error",
        message=message or "no data",
        notes=[],
    )
