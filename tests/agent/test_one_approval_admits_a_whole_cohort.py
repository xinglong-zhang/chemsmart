"""One approval, N authorised elements -- and still one consumption.

The bundle is one-shot, which is right: a human approved one execution of
one plan. But "one-shot" was implemented as "one process", and a wave is
N processes running N members of that same approved plan. Today the
second element either finds the bundle consumed and is admitted through
the *continuation* path -- whose contract says it authorises nothing and
which demands a run state that the startup race may not have written yet
-- or, arriving first, is refused outright as
"a second independent execution of a consumed bundle".

So membership, not a second claim. The cohort manifest already names
exactly which calculations this approval covers and which element runs
which. An element is authorised by *being in it*; the approval is still
consumed once.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.cohort import (
    authorise_cohort_element,
    build_cohort_manifest,
)


def _manifest(tmp_path, nodes=("a1", "a2", "a3")):
    manifest = build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="a" * 64,
        node_ids=nodes,
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    )
    manifest.write(tmp_path)
    return manifest


def test_every_element_of_the_cohort_is_authorised(tmp_path):
    manifest = _manifest(tmp_path)
    for element, node_id in enumerate(manifest.node_ids):
        assert (
            authorise_cohort_element(
                tmp_path, element=element, bundle_sha256="a" * 64
            )
            == node_id
        )


def test_an_element_outside_the_cohort_is_not_authorised(tmp_path):
    _manifest(tmp_path)
    with pytest.raises(ContractError, match="element 7"):
        authorise_cohort_element(tmp_path, element=7, bundle_sha256="a" * 64)


def test_a_different_approval_cannot_borrow_this_cohort(tmp_path):
    """Membership authorises against *this* approval and no other."""

    _manifest(tmp_path)
    with pytest.raises(ContractError, match="approval"):
        authorise_cohort_element(tmp_path, element=0, bundle_sha256="b" * 64)


def test_without_a_manifest_there_is_no_cohort_authority(tmp_path):
    """A single-job run is not a cohort and must keep its own path."""

    assert (
        authorise_cohort_element(
            tmp_path, element=None, bundle_sha256="a" * 64
        )
        is None
    )


def test_the_executor_takes_membership_instead_of_a_second_claim():
    """A claim path nothing changed would still refuse element two."""

    import inspect

    from chemsmart.agent import executor

    source = inspect.getsource(executor)
    assert (
        "authorise_cohort_element" in source or "cohort_element" in source
    ), (
        "the executor still claims the bundle per process, so the second "
        "array element is refused as a second independent execution"
    )
