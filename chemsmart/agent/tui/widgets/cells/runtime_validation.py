"""Runtime validation transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text

from .base import BaseCell

_REMOTE_FIELD_LABELS = {
    "server.queue required": "queue",
    "server.account required": "account",
    "server.scratch_dir required": "scratch_dir",
    "server.modules_or_executable_path required": (
        "modules_or_executable_path"
    ),
}


class RuntimeValidationCell(BaseCell):
    def __init__(self, validation: dict | None) -> None:
        self.validation = validation or {}
        super().__init__(
            self._build_renderable(),
            title="Runtime check",
            classes="runtime-cell",
        )

    def _build_renderable(self):
        local_ok = bool(self.validation.get("local_ok"))
        local_issues = self.validation.get("local_issues") or []
        remote_unknown = self.validation.get("remote_unknown") or []

        if local_ok and remote_unknown and not local_issues:
            count = len(remote_unknown)
            return Group(
                Text("! 원격 제출 정보가 아직 부족합니다", style="bold"),
                Text(
                    f"로컬 OK / 원격 정보 {count}개 필요",
                    style="bold",
                ),
                Text(
                    (
                        "dry-run 결과는 준비됐지만, 실제 제출 전 서버 설정 "
                        "확인이 필요합니다."
                    )
                ),
                Text(_compact_remote_unknown(remote_unknown), style="dim"),
            )

        lines = [Text(_header_text(local_ok, local_issues, remote_unknown))]
        if local_issues:
            lines.append(Text(""))
            lines.append(Text("로컬 환경 확인 필요", style="bold"))
            lines.extend(Text(f"- {issue}") for issue in local_issues)
        if remote_unknown:
            lines.append(Text(""))
            lines.append(Text("확인 필요", style="bold"))
            lines.extend(Text(f"- {issue}") for issue in remote_unknown)
        if not local_issues and not remote_unknown:
            lines.extend(
                [
                    Text(""),
                    Text(
                        (
                            "All runtime checks passed."
                            if self.validation
                            else "No runtime validation data."
                        ),
                        style="dim" if not self.validation else "",
                    ),
                ]
            )
        return Group(*lines)


def _header_text(
    local_ok: bool,
    local_issues: list[str],
    remote_unknown: list[str],
) -> str:
    if local_ok and not local_issues and not remote_unknown:
        return "✓ 로컬/원격 실행 준비됨"
    if local_issues:
        return "✖ 로컬 환경 확인 필요"
    if remote_unknown:
        return f"! 로컬 OK / 원격 정보 {len(remote_unknown)}개 필요"
    return "Runtime 상태를 확인할 수 없습니다"


def _compact_remote_unknown(remote_unknown: list[str]) -> str:
    labels = [
        _REMOTE_FIELD_LABELS.get(item)
        for item in remote_unknown
        if _REMOTE_FIELD_LABELS.get(item)
    ]
    if labels:
        return " / ".join(labels)
    return " / ".join(str(item) for item in remote_unknown)
