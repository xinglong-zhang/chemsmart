"""A human projection of the host-rendered analysis report.

The report file on disk is provenance and stays byte-untouched: completion
receipts, plan digests, and per-claim source receipts all live there. A
chemist reading the terminal needs the science it carries -- conditions,
claims with values and units, and the decision prose -- so this module
re-renders the markdown with the bookkeeping removed. The parser is keyed
to the host renderer's exact headings (a source-pin test fails loudly if
the renderer drifts); unknown sections pass through rather than being
dropped, with bare digests elided from display only.
"""

from __future__ import annotations

import re
from typing import Any

from rich.console import Group
from rich.markdown import Markdown
from rich.table import Table
from rich.text import Text

from chemsmart.agent.report_format import (
    CLAIM_RECORD_LABEL,
    CLAIMS_HEADING,
    COMPLETION_RECEIPT_LABEL,
    CONDITIONS_HEADING,
    DECISION_SECTIONS,
    EVIDENCE_COLUMN,
    HOST_REPORT_TITLE,
    NO_DECISION_PREFIX,
    SOURCE_RECEIPT_COLUMN,
    THERMO_CONDITIONS_HEADING,
    TOOLCHAIN_PLAN_LABEL,
)

_HEX64 = re.compile(r"\b[0-9a-f]{64}\b")
_PROVENANCE_LINE = re.compile(
    rf"^({COMPLETION_RECEIPT_LABEL}|{TOOLCHAIN_PLAN_LABEL}"
    rf"|{CLAIM_RECORD_LABEL}): `[0-9a-f]{{64}}`\s*$"
)
_NO_DECISION_PREFIX = NO_DECISION_PREFIX
_NO_DECISION_HUMAN = (
    "No interpretation was recorded in this run. The numbers above are "
    "host-validated results; the scientific reading remains yours to make "
    "in the session."
)
_DECISION_SECTIONS = set(DECISION_SECTIONS)


def looks_like_host_report(text: str) -> bool:
    return str(text).lstrip().startswith(HOST_REPORT_TITLE)


def _cells(line: str) -> list[str]:
    return [
        cell.strip().strip("`") for cell in line.strip().strip("|").split("|")
    ]


def _is_separator(line: str) -> bool:
    stripped = line.strip()
    return bool(stripped) and set(stripped) <= {"|", "-", ":", " "}


def _read_table(
    lines: list[str], index: int
) -> tuple[list[str], list[list[str]], int]:
    headers = _cells(lines[index])
    index += 1
    if index < len(lines) and _is_separator(lines[index]):
        index += 1
    rows: list[list[str]] = []
    while index < len(lines) and lines[index].strip().startswith("|"):
        rows.append(_cells(lines[index]))
        index += 1
    return headers, rows, index


def _table(
    title: str, headers: list[str], rows: list[list[str]], *, drop: set[int]
) -> Table:
    table = Table(title=title)
    kept = [i for i in range(len(headers)) if i not in drop]
    for position, i in enumerate(kept):
        table.add_column(
            headers[i],
            style="bold cyan" if position == 0 else "",
            justify=(
                "right"
                if headers[i].lower()
                in {"value", "temperature (k)", "pressure (atm)"}
                else "left"
            ),
        )
    for row in rows:
        table.add_row(*[row[i] if i < len(row) else "" for i in kept])
    return table


def render_report_for_humans(markdown_text: str) -> Group:
    lines = str(markdown_text).splitlines()
    blocks: list[Any] = []
    passthrough: list[str] = []

    def flush() -> None:
        if passthrough:
            body = "\n".join(_HEX64.sub("…", line) for line in passthrough)
            if body.strip():
                blocks.append(Markdown(body))
            passthrough.clear()

    index = 0
    while index < len(lines):
        line = lines[index]
        stripped = line.strip()
        if (
            stripped == HOST_REPORT_TITLE.lstrip("# ").strip()
            or stripped == HOST_REPORT_TITLE
        ):
            flush()
            blocks.append(Text("Analysis results", style="bold"))
            index += 1
            continue
        if _PROVENANCE_LINE.match(stripped):
            index += 1
            continue
        if stripped.startswith(_NO_DECISION_PREFIX):
            flush()
            blocks.append(Text(_NO_DECISION_HUMAN, style="italic"))
            # The host wraps this sentence over several lines; consume them.
            index += 1
            while (
                index < len(lines)
                and lines[index].strip()
                and not lines[index].startswith("#")
            ):
                index += 1
            continue
        if stripped == CONDITIONS_HEADING:
            flush()
            index += 1
            while index < len(lines) and not lines[index].strip().startswith(
                "|"
            ):
                index += 1
            headers, rows, index = _read_table(lines, index)
            drop: set[int] = set()
            if EVIDENCE_COLUMN in headers:
                evidence = headers.index(EVIDENCE_COLUMN)
                cells = [row[evidence] for row in rows if evidence < len(row)]
                if cells and all(
                    re.fullmatch(r"[0-9a-f]{64}", cell) for cell in cells
                ):
                    drop = {evidence}
            blocks.append(_table("Conditions", headers, rows, drop=drop))
            continue
        if stripped == THERMO_CONDITIONS_HEADING:
            flush()
            index += 1
            while index < len(lines) and not lines[index].strip().startswith(
                "|"
            ):
                index += 1
            headers, rows, index = _read_table(lines, index)
            blocks.append(
                _table("Thermochemical conditions", headers, rows, drop=set())
            )
            continue
        if stripped == CLAIMS_HEADING:
            flush()
            index += 1
            while index < len(lines) and not lines[index].strip().startswith(
                "|"
            ):
                index += 1
            headers, rows, index = _read_table(lines, index)
            drop = (
                {headers.index(SOURCE_RECEIPT_COLUMN)}
                if SOURCE_RECEIPT_COLUMN in headers
                else set()
            )
            blocks.append(_table("Results", headers, rows, drop=drop))
            continue
        if stripped.startswith("## ") and stripped[3:] in _DECISION_SECTIONS:
            flush()
            section = [stripped[3:]]
            index += 1
            bullets = []
            while index < len(lines) and not lines[index].startswith("## "):
                if lines[index].strip():
                    bullets.append(lines[index])
                index += 1
            blocks.append(
                Markdown(f"**{section[0]}**\n\n" + "\n".join(bullets))
            )
            continue
        passthrough.append(line)
        index += 1
    flush()
    return Group(*blocks)


__all__ = [
    "HOST_REPORT_TITLE",
    "looks_like_host_report",
    "render_report_for_humans",
]
