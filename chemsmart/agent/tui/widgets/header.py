"""Chemsmart chat header wordmark widget."""

from __future__ import annotations

import re

from rich.style import Style
from rich.text import Text
from textual.color import Color, Gradient
from textual.widgets import Static


class ChemsmartHeader(Static):
    """Render a compact chemsmart wordmark for the chat screen."""

    DEFAULT_CSS = """
    ChemsmartHeader {
        height: 1;
    }
    """

    plain_wordmark = "CHEMSMART"
    styled_wordmark = "◢ CHEM▸SMART"

    def render(self) -> Text:
        if getattr(self.app, "plain", False):
            return Text(self.plain_wordmark, end="")
        return self._render_styled_wordmark()

    def _render_styled_wordmark(self) -> Text:
        primary = self._theme_color("primary", "#0178D4")
        accent = self._theme_color("accent", "#FFA62B")
        smart_gradient = Gradient.from_colors(primary, accent, quality=24)

        text = Text(no_wrap=True, overflow="ellipsis", end="")
        text.append("◢", Style(color=primary.rich_color, bold=True))
        text.append(" ", Style(dim=True))
        text.append("CHEM", Style(color=primary.rich_color, bold=True))
        text.append("▸", Style(color=accent.rich_color, bold=True))
        self._append_gradient_text(text, "SMART", smart_gradient)
        return text

    def _append_gradient_text(
        self, text: Text, value: str, gradient: Gradient
    ) -> None:
        last_index = max(len(value) - 1, 1)
        for index, character in enumerate(value):
            position = index / last_index
            text.append(
                character,
                Style(color=gradient.get_rich_color(position), bold=True),
            )

    def _theme_color(self, name: str, fallback: str) -> Color:
        if self.app is None:
            return Color.parse(fallback)
        variables = self.app.get_css_variables()
        color_value = variables.get(name, fallback)
        return Color.parse(color_value)

    @classmethod
    def normalize_wordmark(cls, value: str) -> str:
        """Normalize the rendered wordmark for tests."""

        return re.sub(r"[^A-Za-z]", "", value).upper()
