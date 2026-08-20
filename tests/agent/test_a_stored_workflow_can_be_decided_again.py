"""Deciding twice on the same science must be possible, and must be a decision.

A review is inert and re-loadable, and neither its digest nor a bundle's carries
a run directory, a session id, or a clock -- so the same stored workflow
re-presents identically on an undrifted host.  Re-running approved work was
still effectively blocked, because the default approval id is derived from the
review digest alone: a second honest decision reproduced the first one's
identity, its resolution event conflicted with the persisted one, and the run
aborted before any bundle was written.  A campaign cycle recorded its owed
re-observation as blocked for exactly this reason, with two decision records
that differ only in their timestamps.

These tests execute that path.  The one-shot rule is untouched throughout: a
spent approval stays spent, and what becomes possible is taking a *new*
decision, not reusing the old one.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.execution import workflow_execution_review_json
from chemsmart.agent.live_session import (
    claim_workflow_execution_approval_bundle,
    inspect_workflow_execution_replay,
    replay_approval_id,
    resolve_workflow_execution_review,
)
from tests.agent.test_exact_execution_approval_chain import _review


def _stored_review(tmp_path):
    review = _review(tmp_path)
    path = tmp_path / "review.json"
    path.write_text(workflow_execution_review_json(review), encoding="utf-8")
    return review, path


def _approve(review, review_path, *, approval_id, scope):
    scope.mkdir(parents=True, exist_ok=True)
    return resolve_workflow_execution_review(
        review_file=review_path,
        reviewed_sha256=review.review_sha256,
        decision="approve",
        actor="human",
        output_file=scope / "bundle.json",
        decision_log=scope / "decisions.jsonl",
        approval_id=approval_id,
    )


def test_a_second_decision_needs_its_own_identity_and_its_own_log(tmp_path):
    """The blocker, reproduced and then cleared."""

    review, review_path = _stored_review(tmp_path)
    workspace = review.request.workspace

    _, first = _approve(
        review, review_path, approval_id="approval-first", scope=tmp_path / "a"
    )
    claim_workflow_execution_approval_bundle(first, workspace=workspace)

    # Re-deciding into the same scope collides on the resolution event, which
    # is what the default id did, and no bundle is produced at all.
    with pytest.raises(ContractError):
        _approve(
            review,
            review_path,
            approval_id="approval-second",
            scope=tmp_path / "a",
        )

    # Its own identity and its own log: a real second decision.
    _, second = _approve(
        review,
        review_path,
        approval_id="approval-second",
        scope=tmp_path / "b",
    )
    assert second is not None
    assert second.bundle_sha256 != first.bundle_sha256
    claim_workflow_execution_approval_bundle(second, workspace=workspace)


def test_the_spent_approval_stays_spent(tmp_path):
    """Nothing here relaxes one-shot consumption."""

    review, review_path = _stored_review(tmp_path)
    workspace = review.request.workspace
    _, bundle = _approve(
        review, review_path, approval_id="approval-once", scope=tmp_path / "a"
    )

    claim_workflow_execution_approval_bundle(bundle, workspace=workspace)
    with pytest.raises(ContractError, match="already been consumed"):
        claim_workflow_execution_approval_bundle(bundle, workspace=workspace)


def test_the_preview_names_the_identity_already_burned(tmp_path):
    """ "Already consumed" becomes an instruction instead of a dead end."""

    review, review_path = _stored_review(tmp_path)
    workspace = review.request.workspace
    _, bundle = _approve(
        review, review_path, approval_id="approval-once", scope=tmp_path / "a"
    )
    claim_workflow_execution_approval_bundle(bundle, workspace=workspace)

    report = inspect_workflow_execution_replay(
        review_file=review_path, workspace=workspace
    )

    assert report["review_sha256"] == review.review_sha256
    assert "approval-once" in report["previously_consumed_approval_ids"]


def test_a_replay_identity_does_not_collide_with_the_derived_default():
    first, second = replay_approval_id(), replay_approval_id()
    assert first != second
    assert first.startswith("replay-")


def test_absent_approved_bytes_are_named_so_no_decision_is_offered(tmp_path):
    """Approving a certain dead end would burn a decision for nothing.

    The fixture's workspace holds no artifact bytes at all, which is the same
    condition as an approved input having been moved away: every binding the
    review names is unresolvable, so the report says so and the command refuses
    to offer a decision on it.
    """

    review, review_path = _stored_review(tmp_path)

    report = inspect_workflow_execution_replay(
        review_file=review_path, workspace=review.request.workspace
    )

    assert report["missing_approved_artifacts"]
    assert not report["approved_artifacts_present"]


def test_a_review_for_another_workspace_is_refused(tmp_path):
    review, review_path = _stored_review(tmp_path)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()

    with pytest.raises(ContractError, match="another workspace"):
        inspect_workflow_execution_replay(
            review_file=review_path, workspace=elsewhere
        )
