"""Receipts the host produced must not be handed back to it by the model.

`synthesize_command` and `preflight_program_node` each take optional arguments
naming receipts the host itself created moments earlier, for the same
invocation, and each shipped undescribed:

* `synthesize_command.project_validation_receipt_sha256`
* `preflight_program_node.project_validation_receipt_sha256`
* `preflight_program_node.safe_preview_receipt_sha256`

Omitting them did not produce a message about the omission. Omitting the first
made the compiler say

    project 'X' is bound but has no validation receipt; call
    validate_project_yaml on it first

to a caller that had just called it successfully. Omitting the others made the
preflight `blocked`, and execution then refused three layers away with "node
requires a green safe-preview preflight", naming neither.

Measured by driving an approved plan through the real host: four of the six
blockers on the approval-to-execution path were of this kind.

The division of labour is the point. The model decides which project and which
node; which receipt records that decision is bookkeeping with exactly one right
answer, and the host holds it. An explicitly supplied digest still wins, so a
caller that wants to name a specific receipt can.
"""

import inspect

from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def test_synthesize_resolves_the_validation_receipt_it_already_has():
    source = inspect.getsource(CommandCompiledToolHostV1._synthesize_command)
    assert "_resolve_project_validation" in source
    assert (
        "if validation_digest" in source
    ), "an explicitly supplied receipt must still take precedence"


def test_preflight_resolves_the_receipts_the_host_just_produced():
    source = inspect.getsource(
        CommandCompiledToolHostV1._preflight_program_node
    )
    assert "invocation.project_receipt_sha256" in source, (
        "the invocation records the validation receipt it was compiled "
        "against, so there is nothing to re-derive"
    )
    assert "_resolve_safe_preview" in source
    assert "if project_digest" in source and "if preview_digest" in source


def test_the_safe_preview_lookup_is_keyed_by_invocation():
    source = inspect.getsource(CommandCompiledToolHostV1._resolve_safe_preview)
    assert "invocation_sha256 == invocation_sha256" in source.replace(
        "receipt.", ""
    )


def test_an_ambiguous_validation_is_refused_rather_than_guessed():
    """Silently picking one of several valid receipts would hide a real fork."""
    source = inspect.getsource(
        CommandCompiledToolHostV1._resolve_project_validation
    )
    assert "len(matches) > 1" in source
    assert "pass project_validation_receipt_sha256 to say" in source


def test_resolution_returns_none_so_the_real_missing_case_still_speaks():
    source = inspect.getsource(
        CommandCompiledToolHostV1._resolve_project_validation
    )
    assert "return matches[0] if matches else None" in source, (
        "when no receipt exists the caller must still get the message that "
        "validation is genuinely missing"
    )


def test_the_predicate_matches_the_one_prepare_program_node_uses():
    """Two resolvers disagreeing would be worse than one being absent."""
    resolver = inspect.getsource(
        CommandCompiledToolHostV1._resolve_project_validation
    )
    for clause in (
        "receipt.project_artifact_id == project.artifact_id",
        "receipt.project_sha256 == project.sha256",
        "receipt.capability_receipt_sha256 == capability.receipt_sha256",
        'receipt.status == "valid"',
    ):
        assert clause in resolver, clause
