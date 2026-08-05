"""A rejection must carry what the caller needs to retry differently.

Two live runs stalled on this. The first repeated a DAG rejection verbatim
because it named neither the node nor the value; the second collided five
times on "artifact ID is already registered" because it never said which IDs
were taken. A caller that cannot act on a rejection retries unchanged, which
looks like a model defect and is a message defect.
"""

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def _get(values, key, label):
    return CommandCompiledToolHostV1._get(values, key, label)


def test_an_unknown_id_names_what_was_asked_for_and_what_exists():
    with pytest.raises(ContractError) as excinfo:
        _get({"alpha": 1, "beta": 2}, "gamma", "project render receipt")
    message = str(excinfo.value)
    assert "gamma" in message, "the rejected ID must be named"
    assert "alpha" in message and "beta" in message, "bound IDs must be listed"


def test_an_empty_registry_says_so_rather_than_listing_nothing():
    with pytest.raises(ContractError, match="no run receipt is bound yet"):
        _get({}, "anything", "run receipt")


def test_a_large_registry_is_truncated_but_still_counted():
    values = {f"id{index:02d}": index for index in range(20)}
    with pytest.raises(ContractError) as excinfo:
        _get(values, "missing", "settings object")
    message = str(excinfo.value)
    assert "id00" in message
    assert "12 more" in message, "the caller should learn the list is partial"


def test_a_duplicate_project_id_lists_the_taken_ids():
    """The exact rejection a live run hit five times in a row."""

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    host.approved_workspace = object()
    host.artifacts = {"project-a": object(), "project-b": object()}
    with pytest.raises(ContractError) as excinfo:
        host._promote_project_yaml(
            "turn", {"artifact_id": "project-a", "render_receipt_sha256": "x"}
        )
    message = str(excinfo.value)
    assert "project-a" in message
    assert "project-b" in message, "the caller must see what is already taken"
