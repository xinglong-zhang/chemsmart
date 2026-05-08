"""Dry-run input transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text

from .base import BaseCell


class DryRunInputCell(BaseCell):
    def __init__(self, content: str, inputfile: str | None = None) -> None:
        header = Text()
        if inputfile:
            header.append("Input file: ", style="bold")
            header.append(inputfile)
        else:
            header.append("Dry-run input generated", style="bold")
        body = Text(content)
        super().__init__(
            Group(header, Text(""), body),
            title="Dry run input",
            classes="dry-run-cell",
        )
