"""Explicit workspace project-YAML write confirmation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.screen import ModalScreen
from textual.widgets import Button, Static


@dataclass(slots=True, frozen=True)
class ProjectWriteResult:
    action: str
    project_name: str


class ProjectWriteOverlay(ModalScreen[ProjectWriteResult | None]):
    """Ask whether a validated candidate replaces the active YAML."""

    BINDINGS = [
        Binding("y", "overwrite", "Overwrite", show=False, priority=True),
        Binding("n", "write_new", "Write new", show=False, priority=True),
        Binding("escape", "cancel", "Cancel", show=False, priority=True),
    ]

    DEFAULT_CSS = """
    ProjectWriteOverlay {
        align: center bottom;
    }

    #project-write-modal {
        width: 88;
        max-width: 96%;
        height: auto;
        margin-bottom: 7;
        border: round $warning;
        padding: 1 2;
        background: $surface;
    }

    #project-write-actions {
        width: 100%;
        height: auto;
        margin-top: 1;
        align-horizontal: center;
    }

    #project-write-actions Button {
        margin: 0 1;
    }
    """

    def __init__(
        self,
        *,
        target: Path,
        overwrite_project: str,
        new_project: str,
        target_exists: bool,
    ) -> None:
        super().__init__()
        self.target = target
        self.overwrite_project = overwrite_project
        self.new_project = new_project
        self.target_exists = target_exists

    def compose(self) -> ComposeResult:
        if self.target_exists:
            question = (
                f"Shell: overwrite [{self.target.name}]?\n\n"
                "The candidate has already passed the project-YAML validator. "
                "Choose whether it replaces the active settings or is saved "
                "as a separate workspace project."
            )
            overwrite_label = "Yes, overwrite settings (Y)"
            new_label = f"No, write {self.new_project}.yaml (N)"
        else:
            question = (
                f"Shell: write [{self.target.name}]?\n\n"
                "The candidate has passed the project-YAML validator and will "
                "become the active workspace project."
            )
            overwrite_label = "Yes, write settings (Y)"
            new_label = "No, cancel (N)"
        with Vertical(id="project-write-modal"):
            yield Static(question, id="project-write-summary")
            with Horizontal(id="project-write-actions"):
                yield Button(
                    overwrite_label,
                    id="project-write-overwrite",
                    variant="warning" if self.target_exists else "success",
                )
                yield Button(new_label, id="project-write-new")

    def on_mount(self) -> None:
        self.query_one("#project-write-summary", Static).border_title = (
            "Workspace project YAML"
        )

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "project-write-overwrite":
            self.action_overwrite()
        elif event.button.id == "project-write-new":
            self.action_write_new()

    def action_overwrite(self) -> None:
        self.dismiss(
            ProjectWriteResult(
                "overwrite" if self.target_exists else "write",
                self.overwrite_project,
            )
        )

    def action_write_new(self) -> None:
        if not self.target_exists:
            self.dismiss(None)
            return
        self.dismiss(ProjectWriteResult("write_new", self.new_project))

    def action_cancel(self) -> None:
        self.dismiss(None)
