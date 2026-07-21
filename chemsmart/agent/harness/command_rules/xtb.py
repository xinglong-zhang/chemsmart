"""Deterministic contract rules for xTB commands."""

from __future__ import annotations

from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    reject,
)
from chemsmart.agent.harness.command_rules.tokens import has_option

SOLVENT_MODEL_ALIASES = ("-sm", "--solvent-model")
SOLVENT_ID_ALIASES = ("-si", "--solvent-id")


def solvent_contract_issues(
    program_tokens: list[str],
) -> tuple[CommandContractIssue, ...]:
    """Reject half-specified xTB solvation.

    xTB needs both the model (``alpb``, ``gbsa``, ...) and the solvent
    identifier; either one alone silently runs in the gas phase, which is a
    different calculation from the one that was asked for. Click cannot
    express the pairing, so the contract owns it.
    """

    has_model = has_option(program_tokens, SOLVENT_MODEL_ALIASES)
    has_solvent = has_option(program_tokens, SOLVENT_ID_ALIASES)
    if has_model == has_solvent:
        return ()
    missing = (
        "xTB solvent identifier (-si)"
        if has_model
        else "xTB solvent model (-sm)"
    )
    return (
        reject(
            "cmd.contract.xtb_solvent_pair",
            (
                "xTB solvation requires both --solvent-model and "
                "--solvent-id; specifying one alone runs in the gas phase"
            ),
            {
                "program": "xtb",
                "solvent_model_present": has_model,
                "solvent_id_present": has_solvent,
            },
            (missing,),
        ),
    )


__all__ = [
    "SOLVENT_ID_ALIASES",
    "SOLVENT_MODEL_ALIASES",
    "solvent_contract_issues",
]
