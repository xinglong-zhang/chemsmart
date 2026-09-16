"""A contract nothing reads is a paragraph.

Two of this round's own functions shipped with no production reader, and
the suite was green both times because each had tests that called it
directly. `cohort_completion` -- the host's answer to "has this wave
ended" -- decided nothing, while a file's existence decided the wake;
`execution_result_file` owned the per-element path while the array
script spelled its own.

Being *called by a test* is the disguise. This asks the only question
that separates a contract from a paragraph: does anything outside this
module, in the shipped package, name it?
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

_MODULE = Path("chemsmart/agent/cohort.py")


def _public_functions(path: Path) -> tuple[str, ...]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return tuple(
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and not node.name.startswith("_")
    )


def _package_readers(name: str) -> tuple[str, ...]:
    """Every shipped module that names this symbol, excluding its home."""

    found = []
    for path in Path("chemsmart").rglob("*.py"):
        if path == _MODULE:
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (OSError, SyntaxError):  # pragma: no cover
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and node.id == name:
                found.append(str(path))
                break
            if isinstance(node, ast.Attribute) and node.attr == name:
                found.append(str(path))
                break
        else:
            continue
    return tuple(sorted(set(found)))


@pytest.mark.parametrize("name", _public_functions(_MODULE))
def test_every_cohort_function_decides_something(name):
    readers = _package_readers(name)
    assert readers, (
        f"{name} is declared in {_MODULE} and read by nothing the "
        "package ships, so whatever it decides is decided somewhere "
        "else -- which is how the wake came to be gated on a filename "
        "while the barrier sat unused beside it"
    )


def test_the_check_can_fail():
    """A guard that cannot go red guards nothing."""

    assert _package_readers("a_name_this_package_does_not_contain") == ()
