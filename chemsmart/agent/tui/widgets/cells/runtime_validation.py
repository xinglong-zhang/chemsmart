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
                Text("! Remote submission info is incomplete", style="bold"),
                Text(
                    f"Local OK / {count} remote item(s) needed",
                    style="bold",
                ),
                Text(
                    (
                        "Dry-run results are ready, but server settings need "
                        "to be confirmed before actual submission."
                    )
                ),
                Text(_compact_remote_unknown(remote_unknown), style="dim"),
            )

        lines = [Text(_header_text(local_ok, local_issues, remote_unknown))]
        if local_issues:
            lines.append(Text(""))
            lines.append(Text("Local environment needs review", style="bold"))
            lines.extend(Text(f"- {issue}") for issue in local_issues)
        if remote_unknown:
            lines.append(Text(""))
            lines.append(Text("Needs review", style="bold"))
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
        return "✓ Local/remote run ready"
    if local_issues:
        return "✖ Local environment needs review"
    if remote_unknown:
        return f"! Local OK / {len(remote_unknown)} remote item(s) needed"
    return "Runtime status unavailable"


def _compact_remote_unknown(remote_unknown: list[str]) -> str:
    labels = [
        _REMOTE_FIELD_LABELS.get(item)
        for item in remote_unknown
        if _REMOTE_FIELD_LABELS.get(item)
    ]
    if labels:
        return " / ".join(labels)
    return " / ".join(str(item) for item in remote_unknown)
