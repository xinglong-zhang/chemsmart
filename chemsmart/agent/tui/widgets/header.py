"""Chemsmart chat header wordmark widget."""

from __future__ import annotations

import re

from rich.style import Style
from rich.text import Text
from textual.color import Color, ColorParseError
from textual.widgets import Static


class ChemsmartHeader(Static):
    """Render a compact chemsmart wordmark for the chat screen."""

    DEFAULT_CSS = """
    ChemsmartHeader {
        height: 1;
    }
    """

    plain_wordmark = "CHEMSMART"
    styled_wordmark = "chem · smart"

    def render(self) -> Text:
        if getattr(self.app, "plain", False):
            return Text(self.plain_wordmark, end="")
        return self._render_styled_wordmark()

    def _render_styled_wordmark(self) -> Text:
        body_style = Style(
            color=self._theme_color(
                "$text-muted", fallback="#B8B8B8"
            ).rich_color,
        )
        divider_style = Style(
            color=self._theme_color(
                "$foreground",
                "$text-disabled",
                fallback="#E0E0E0",
            ).rich_color,
        )

        text = Text(
            self.styled_wordmark,
            style=body_style,
            no_wrap=True,
            overflow="ellipsis",
            end="",
        )
        divider_index = self.styled_wordmark.index("·")
        text.stylize(divider_style, divider_index, divider_index + 1)
        return text

    def _theme_color(self, *names: str, fallback: str) -> Color:
        if self.app is None:
            return Color.parse(fallback)

        variables = self.app.get_css_variables()
        for name in names:
            for candidate in self._theme_candidates(name):
                color_value = variables.get(candidate)
                if color_value is None:
                    continue
                try:
                    return Color.parse(color_value)
                except ColorParseError:
                    continue
        return Color.parse(fallback)

    def _theme_candidates(self, name: str) -> tuple[str, ...]:
        token = name.removeprefix("$")
        if token == "text":
            return (token, "foreground")
        if token.startswith("text-"):
            return (token, f"foreground-{token.removeprefix('text-')}")
        return (token,)

    @classmethod
    def normalize_wordmark(cls, value: str) -> str:
        """Normalize the rendered wordmark for tests."""

        return re.sub(r"[^A-Za-z]", "", value).upper()
