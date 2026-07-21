"""Invariant checks over the argv xTB was actually invoked with."""

from __future__ import annotations

from typing import Any

from chemsmart.agent.harness.generated_common import reject
from chemsmart.agent.harness.models import InvariantIssue

#: ``--gfnff`` is a standalone flag; the numbered Hamiltonians are selected as
#: ``--gfn <n>``. Source: ``chemsmart/jobs/xtb/runner.py``.
_GFN_FF = "gfnff"


def solvent_issues(
    route: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    """Reject a rendered call that drops or half-applies requested solvation."""

    model = chemistry.get("solvent_model")
    solvent = chemistry.get("solvent_id")
    if not model and not solvent:
        return []
    if not model or not solvent:
        return [
            reject(
                "input.xtb.solvent_pair",
                (
                    "xTB solvation needs both a solvent model and a solvent "
                    "identifier"
                ),
                {
                    **evidence,
                    "requested_solvent_model": model,
                    "requested_solvent_id": solvent,
                },
            )
        ]
    tokens = route.split()
    flag = f"--{str(model).lower()}"
    if flag in tokens and _value_after(tokens, flag) == str(solvent).lower():
        return []
    return [
        reject(
            "input.xtb.solvent_preservation",
            (
                "generated xTB command does not preserve the requested "
                "implicit solvation"
            ),
            {
                **evidence,
                "requested_solvent_model": model,
                "requested_solvent_id": solvent,
                "expected_flag": flag,
            },
        )
    ]


def gfn_version_issues(
    route: str,
    chemistry: dict[str, Any],
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    """Reject a rendered call that silently changes the GFN Hamiltonian."""

    requested = chemistry.get("gfn_version")
    if not requested:
        return []
    requested = str(requested).lower()
    tokens = route.split()
    if requested == _GFN_FF:
        rendered_ok = "--gfnff" in tokens
    else:
        rendered_ok = (
            "--gfn" in tokens
            and _value_after(tokens, "--gfn") == requested[-1:]
        )
    if rendered_ok:
        return []
    return [
        reject(
            "input.xtb.gfn_preservation",
            "generated xTB command does not preserve the requested GFN method",
            {**evidence, "requested_gfn_version": requested},
        )
    ]


def _value_after(tokens: list[str], flag: str) -> str | None:
    if flag not in tokens:
        return None
    index = tokens.index(flag)
    if index + 1 >= len(tokens):
        return None
    return tokens[index + 1].lower()


__all__ = ["gfn_version_issues", "solvent_issues"]
