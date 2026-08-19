"""A producer-fed analysis node could never be reported complete.

The receipt matcher routed every dependency through its analysis-receipt
ledger, and calculation nodes are never in that ledger -- so an extraction
node consuming a calculation's output stayed "waiting" forever, even after
the producer had executed and validated and the extraction receipt existed.
This is the session-and-executor truth source for analysis completion, so
the whole approved-analysis-chain feature would have reported nothing done.
Calculation dependencies are now satisfied by validated execution receipts,
and a program_output input matches extraction receipts over the registered
result artifacts of its producer.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    build_scientific_toolchain_plan,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)

_RESULT = Path("tests/data/ORCATests/outputs/CO2.out")


def _plan():
    calculation = CommandNodeIntentV1(
        node_id="sp",
        program="orca",
        jobtype="sp",
        project_role="r",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geometry",
                artifact_class="geometry_xyz",
                artifact_id="start",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id="sp-out", artifact_class="orca_output"
            ),
        ),
        unresolved_fields=(),
    )
    analysis = AnalysisNodeIntentV1(
        node_id="extract-sp",
        analysis_kind="result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )
    return build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(calculation,),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(analysis,),
        required_output_ids=("e",),
    )


def _host_with_extraction(tmp_path):
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    resolved = _RESULT.resolve()
    artifact = TrustedArtifactRefV1(
        artifact_id="result.sp.1",
        kind="orca_output",
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )
    host.artifacts[artifact.artifact_id] = artifact
    host.dispatch(
        turn_id="t1",
        tool_name="extract_result_quantities",
        arguments={
            "program": "orca",
            "artifact_id": "result.sp.1",
            "selectors": [{"quantity_id": "e", "selector": "energy"}],
        },
    )
    return host


def test_a_validated_producer_completes_its_extraction_node(tmp_path):
    host = _host_with_extraction(tmp_path)
    host.execution_receipts["sp"] = SimpleNamespace(validated=True)

    matched = host._scientific_toolchain_analysis_receipts(
        _plan(), task_spec_sha256="a" * 64
    )

    assert "extract-sp" in matched
    assert matched["extract-sp"]


def test_an_unvalidated_producer_completes_nothing(tmp_path):
    host = _host_with_extraction(tmp_path)

    matched = host._scientific_toolchain_analysis_receipts(
        _plan(), task_spec_sha256="a" * 64
    )

    assert "extract-sp" not in matched
