"""Focused checks for agent-bound process resource evidence."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _process_observation_findings,
    _write_host_execution_artifact,
)
from chemsmart.utils.process_observation import (
    launch_failure_observation,
    observe_process,
)


def _completed_observation():
    process = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(0.1)"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    return observe_process(
        process,
        timeout_seconds=600,
        memory_limit_mb=4096,
        sample_interval_seconds=0.02,
    ).observation


def test_live_evaluation_binds_resource_observation_without_false_finding():
    observation = _completed_observation()

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="sp",
        charge=0,
        multiplicity=1,
        output_artifacts=(),
        exit_status=0,
        process_observation=observation,
    )

    assert evaluation.observations["process_observation"][
        "receipt_sha256"
    ] == observation.receipt_sha256
    assert not any(
        finding.startswith("execution.process.")
        for finding in evaluation.findings
    )


def test_launch_failure_is_a_deterministic_execution_finding():
    observation = launch_failure_observation(
        timeout_seconds=600,
        memory_limit_mb=4096,
        error_type="FileNotFoundError",
    )

    assert _process_observation_findings(observation) == (
        "execution.process.launch_failed",
    )


def test_resource_receipt_is_an_immutable_execution_output(tmp_path):
    observation = _completed_observation()
    path = tmp_path / "execution-resource.receipt.json"
    payload = json.dumps(
        observation.as_dict(), sort_keys=True, separators=(",", ":")
    ) + "\n"
    _write_host_execution_artifact(path, payload)
    host = object.__new__(CommandCompiledToolHostV1)
    host.artifacts = {}

    artifacts = host._execution_output_artifacts("sp-initial", tmp_path)

    resource = next(
        item
        for item in artifacts
        if Path(item.path).name == "execution-resource.receipt.json"
    )
    assert resource.kind == "json"
    assert resource.sha256
    with pytest.raises(ContractError, match="reserved host artifact"):
        _write_host_execution_artifact(path, payload)
