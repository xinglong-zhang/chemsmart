"""`--cohort-element` without a cohort must not admit anything.

The executor takes cohort membership *instead of* the workspace-level
one-shot claim -- which is right, because the manifest names exactly
which calculations this approval covers and is digest-bound to it. But
the branch ran on the flag alone:

    if cohort_element is not None:
        authorise_cohort_element(...)   # return discarded
        self._bundle_claimed = True
        return

and `authorise_cohort_element` returns `None`, without raising, when
there is no manifest. So `chemsmart agent run --approval-file <bundle>
--run-directory <a fresh directory> --cohort-element 0` skips the
one-shot claim entirely, and `_cohort_scope()` -- also `None` for a
missing manifest -- leaves the walk unbounded. One human approval, a new
run directory each time, and the whole partition runs again.

It is a public, documented flag, and its help text says "Omit it to run
the whole approved partition" -- which is what giving it did. Unlike the
concurrency defects this one is reachable by a person typing a plausible
command.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.cohort import build_cohort_manifest


def _executor(tmp_path, *, element, manifest=True):
    from chemsmart.agent.executor import ApprovedWorkflowExecutor

    run_directory = tmp_path / "run"
    run_directory.mkdir(parents=True, exist_ok=True)
    if manifest:
        build_cohort_manifest(
            goal_id="g1",
            cycle=1,
            bundle_sha256="e" * 64,
            node_ids=("a1", "a2"),
            max_concurrent_tasks=4,
            created_at="2026-09-16T00:00:00+00:00",
        ).write(run_directory)
    return ApprovedWorkflowExecutor(
        host=SimpleNamespace(
            verify_reviewed_real_execution_argv=lambda **_kw: None
        ),
        plan=SimpleNamespace(
            workflow_id="w",
            plan_sha256="b" * 64,
            nodes=(
                SimpleNamespace(node_id="a1", program="orca"),
                SimpleNamespace(node_id="a2", program="orca"),
            ),
        ),
        approval=SimpleNamespace(node_bindings=()),
        frozen_approval=SimpleNamespace(approval_sha256="c" * 64),
        initial_artifacts={},
        project_artifacts=(),
        task_spec_sha256="a" * 64,
        run_directory=run_directory,
        execution_bundle=SimpleNamespace(
            non_executable_node_ids=(),
            bundle_sha256="e" * 64,
            node_review=lambda node_id: SimpleNamespace(node_id=node_id),
        ),
        approval_workspace=tmp_path / "workspace",
        claim_workspace_bundle=False,
        cohort_element=element,
    )


def test_an_element_without_a_cohort_is_refused(tmp_path):
    executor = _executor(tmp_path, element=0, manifest=False)
    with pytest.raises(ContractError, match="cohort"):
        executor._verify_launch_and_claim_once(
            node_id="a1", invocation_sha256="9" * 64
        )
    assert not getattr(executor, "_bundle_claimed", False), (
        "the one-shot claim was marked taken by a run that took no "
        "claim at all, so one approval buys unbounded re-execution"
    )


def test_an_element_of_a_real_cohort_is_admitted(tmp_path):
    executor = _executor(tmp_path, element=1)
    executor._verify_launch_and_claim_once(
            node_id="a1", invocation_sha256="9" * 64
        )
    assert executor._bundle_claimed is True
    assert executor._cohort_scope() == ("a2",)


def test_an_element_outside_the_cohort_is_refused(tmp_path):
    executor = _executor(tmp_path, element=7)
    with pytest.raises(ContractError):
        executor._verify_launch_and_claim_once(
            node_id="a1", invocation_sha256="9" * 64
        )
