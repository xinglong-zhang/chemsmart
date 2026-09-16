"""A cohort is a named, immutable set of approved calculations.

An array index is scheduler representation. It must never become
scientific identity, and nothing durable may key a calculation by it: the
manifest is the only mapping from element to node, so a reader asking
"what did element 3 compute" gets an approved node id and not a position.

Membership is fixed when the cohort is dispatched. A wave is what the
Agent chose to see before reasoning again, so a set that could grow after
submission would be a different experiment than the one the barrier is
waiting for.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.cohort import (
    CohortManifestV1,
    build_cohort_manifest,
    read_cohort_manifest,
)


def _manifest(**overrides):
    values = dict(
        goal_id="g1",
        cycle=1,
        bundle_sha256="a" * 64,
        node_ids=("conformer-1-opt", "conformer-2-opt", "conformer-3-opt"),
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    )
    values.update(overrides)
    return build_cohort_manifest(**values)


def test_an_element_resolves_to_an_approved_node_not_a_position():
    manifest = _manifest()
    assert manifest.node_for_element(0) == "conformer-1-opt"
    assert manifest.node_for_element(2) == "conformer-3-opt"
    assert manifest.element_for_node("conformer-2-opt") == 1


def test_an_element_outside_the_cohort_is_refused_by_name():
    manifest = _manifest()
    with pytest.raises(ContractError, match="element 5"):
        manifest.node_for_element(5)
    with pytest.raises(ContractError, match="not a member"):
        manifest.element_for_node("conformer-9-opt")


def test_the_cohort_id_is_its_content():
    """Two cohorts of the same members in the same order are the same
    cohort; a different order or a different member is a different one."""

    first = _manifest()
    assert first.cohort_id == _manifest().cohort_id
    reordered = _manifest(
        node_ids=("conformer-2-opt", "conformer-1-opt", "conformer-3-opt")
    )
    assert reordered.cohort_id != first.cohort_id
    smaller = _manifest(node_ids=("conformer-1-opt", "conformer-2-opt"))
    assert smaller.cohort_id != first.cohort_id


def test_membership_cannot_be_edited_after_it_is_written(tmp_path):
    manifest = _manifest()
    path = manifest.write(tmp_path)
    again = read_cohort_manifest(tmp_path)
    assert again == manifest

    tampered = path.read_text(encoding="utf-8").replace(
        "conformer-3-opt", "conformer-9-opt"
    )
    path.write_text(tampered, encoding="utf-8")
    with pytest.raises(ContractError, match="digest"):
        read_cohort_manifest(tmp_path)


def test_a_cohort_holds_no_duplicate_and_no_empty_member():
    with pytest.raises(ContractError, match="once"):
        _manifest(node_ids=("a-opt", "a-opt"))
    with pytest.raises(ContractError):
        _manifest(node_ids=())


def test_a_missing_manifest_is_absence_not_an_error(tmp_path):
    """A single-job dispatch writes none, and that is not a failure."""

    assert read_cohort_manifest(tmp_path) is None


def test_the_manifest_records_what_bounds_it():
    manifest = _manifest(max_concurrent_tasks=4)
    record = manifest.public_record()
    assert record["max_concurrent_tasks"] == 4
    assert record["node_ids"] == list(manifest.node_ids)
    assert record["cohort_sha256"] == manifest.cohort_id
    assert isinstance(manifest, CohortManifestV1)
