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

#: The modules this round added or rebuilt. `cohort.py` is where the
#: defect was found twice; `scheduler_job.py` is where it was found a
#: third time, in `parse_squeue_array` and `squeue_array_command`, which
#: were written for the barrier's own question and read by nothing until
#: the poller was wired to them.
_MODULES = (
    Path("chemsmart/agent/cohort.py"),
    Path("chemsmart/settings/probe/scheduler_job.py"),
    Path("chemsmart/settings/scheduler_request.py"),
)


def _public_functions(path: Path) -> tuple[str, ...]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return tuple(
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and not node.name.startswith("_")
    )


def _declared() -> tuple[tuple[str, str], ...]:
    return tuple(
        (str(module), name)
        for module in _MODULES
        for name in _public_functions(module)
    )


def _package_readers(name: str, home: Path | str) -> tuple[str, ...]:
    """Every shipped module that reads this symbol.

    Its own module counts: a helper another function in the same file
    calls is read, and the question is whether *anything* the package
    ships decides with it -- not whether the call crosses a file
    boundary. The definition itself is not a reading of itself, which is
    the only thing excluded.
    """

    found = []
    for path in Path("chemsmart").rglob("*.py"):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (OSError, SyntaxError):  # pragma: no cover
            continue
        for node in ast.walk(tree):
            if (
                isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name == name
            ):
                continue
            if isinstance(node, ast.Name) and node.id == name:
                found.append(str(path))
                break
            if isinstance(node, ast.Attribute) and node.attr == name:
                found.append(str(path))
                break
    del home
    return tuple(sorted(set(found)))


@pytest.mark.parametrize("home,name", _declared())
def test_every_cohort_function_decides_something(home, name):
    readers = _package_readers(name, home)
    assert readers, (
        f"{name} is declared in {home} and read by nothing the "
        "package ships, so whatever it decides is decided somewhere "
        "else -- which is how the wake came to be gated on a filename "
        "while the barrier sat unused beside it"
    )


def test_the_check_can_fail():
    """A guard that cannot go red guards nothing."""

    assert (
        _package_readers(
            "a_name_this_package_does_not_contain", _MODULES[0]
        )
        == ()
    )
