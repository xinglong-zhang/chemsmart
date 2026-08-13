from __future__ import annotations

from pathlib import Path

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.feedback import project_tool_feedback
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.io.native_failure import (
    summarize_gaussian_native_failure,
    summarize_orca_native_failure,
)


_FIXTURES = Path("tests/data/agent/native_failures")


def _lines(name: str) -> tuple[str, ...]:
    return tuple((_FIXTURES / name).read_text(encoding="utf-8").splitlines())


def _artifact(
    name: str, *, artifact_id: str, kind: str
) -> TrustedArtifactRefV1:
    path = (_FIXTURES / name).resolve()
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=kind,
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


def test_orca_parser_classifies_auxiliary_basis_without_returning_input_text():
    summary = summarize_orca_native_failure(_lines("orca_auxiliary_basis.txt"))

    assert summary is not None
    assert summary.termination_state == "error_termination"
    assert summary.error_class == "auxiliary_basis"
    assert len(summary.diagnostic_lines) <= 3
    assert summary.diagnostic_lines == (
        "ORCA rejected the auxiliary basis for a correlated method.",
    )
    assert all("%basis" not in line for line in summary.diagnostic_lines)


def test_orca_parser_uses_bounded_stderr_to_classify_mpi_runtime():
    summary = summarize_orca_native_failure(
        _lines("orca_mpi_output.txt"),
        diagnostic_lines=_lines("orca_mpi_stderr.err"),
    )

    assert summary is not None
    assert summary.error_class == "mpi_runtime"
    assert len(summary.diagnostic_lines) == 1
    assert summary.diagnostic_lines == (
        "ORCA subprocess reported an MPI runtime failure.",
    )


def test_gaussian_parser_keeps_qperr_but_not_route_path_or_secret():
    summary = summarize_gaussian_native_failure(_lines("gaussian_qperr.txt"))

    assert summary is not None
    assert summary.error_class == "input_syntax"
    assert summary.diagnostic_lines == (
        "Gaussian reported QPErr input syntax failure.",
    )
    joined = " ".join(summary.diagnostic_lines)
    assert "def2-TZVP" not in joined
    assert "/opt/gaussian" not in joined
    assert "fixture-secret" not in joined


def test_orca_matched_line_suffix_is_never_exposed_in_public_diagnostics():
    sentinels = (
        "bare-private-token",
        "quoted private value with spaces",
        "bearer-private-value",
        "assigned-private-value",
        "/private/absolute/path",
        "https://private.invalid/diagnostic",
    )
    line = (
        "SCF CONVERGENCE FAILURE "
        + sentinels[0]
        + ' "'
        + sentinels[1]
        + '" Bearer '
        + sentinels[2]
        + " TOKEN="
        + sentinels[3]
        + " "
        + sentinels[4]
        + " "
        + sentinels[5]
    )

    summary = summarize_orca_native_failure(
        (line, "ORCA finished by error termination in SCF")
    )

    assert summary is not None
    assert summary.error_class == "scf_convergence"
    assert summary.diagnostic_lines == ("ORCA SCF did not converge.",)
    public = " ".join(summary.diagnostic_lines)
    assert all(sentinel not in public for sentinel in sentinels)


def test_gaussian_qperr_exposes_only_canonical_template():
    sentinels = (
        "bare-private-token",
        "quoted private value with spaces",
        "bearer-private-value",
        "assigned-private-value",
        "/private/absolute/path",
        "https://private.invalid/diagnostic",
    )
    matched = (
        "QPErr syntax error was detected in the input line "
        + sentinels[0]
        + ' "'
        + sentinels[1]
        + '" Bearer '
        + sentinels[2]
        + " API_KEY="
        + sentinels[3]
        + " "
        + sentinels[4]
        + " "
        + sentinels[5]
    )

    summary = summarize_gaussian_native_failure(
        (matched, "Error termination via Lnk1e in /private/g16/l1.exe")
    )

    assert summary is not None
    assert summary.diagnostic_lines == (
        "Gaussian reported QPErr input syntax failure.",
    )
    public = " ".join(summary.diagnostic_lines)
    assert all(sentinel not in public for sentinel in sentinels)


def test_gaussian_error_link_text_is_never_reflected():
    for untrusted_link in ("LnkAPIKEYSECRET", "Lnk123SECRET"):
        summary = summarize_gaussian_native_failure(
            (
                "QPErr syntax error was detected in the input line",
                f"Error termination via {untrusted_link} in /private/g16.exe",
            )
        )

        assert summary is not None
        assert summary.diagnostic_lines == (
            "Gaussian reported QPErr input syntax failure.",
        )
        assert untrusted_link not in " ".join(summary.diagnostic_lines)


def test_parser_labels_output_without_a_terminal_marker_as_incomplete():
    summary = summarize_orca_native_failure(_lines("orca_truncated.txt"))

    assert summary is not None
    assert summary.termination_state == "incomplete"
    assert summary.error_class == "incomplete_output"
    assert summary.diagnostic_lines == ()


def test_nonzero_orca_execution_is_parsed_and_bound_to_output_artifacts():
    output = _artifact(
        "orca_mpi_output.txt",
        artifact_id="result.orca.output",
        kind="orca_output",
    )
    stderr = _artifact(
        "orca_mpi_stderr.err",
        artifact_id="result.orca.stderr",
        kind="program_output",
    )

    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="opt",
        charge=0,
        multiplicity=1,
        output_artifacts=(output, stderr),
        exit_status=1,
    )

    failure = evaluation.observations["orca"]["native_failure"]
    assert "execution.process.nonzero_or_unknown" in evaluation.findings
    assert "orca.native_failure.mpi_runtime" in evaluation.findings
    assert failure["error_class"] == "mpi_runtime"
    assert [item["artifact_id"] for item in failure["artifact_refs"]] == [
        "result.orca.output",
        "result.orca.stderr",
    ]
    assert all("path" not in item for item in failure["artifact_refs"])


def test_nonzero_gaussian_failure_survives_causal_execution_feedback():
    output = _artifact(
        "gaussian_qperr.txt",
        artifact_id="result.gaussian.output",
        kind="gaussian_output",
    )
    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="gaussian",
        jobtype="opt",
        charge=0,
        multiplicity=1,
        output_artifacts=(output,),
        exit_status=1,
    )
    tool_result = {
        "schema_version": "chemsmart.tool-result.v1",
        "tool": "execute_approved_program_node",
        "status": "ok",
        "result": {
            "execution": {
                "execution_state": "failed",
                "validated": False,
                "findings": evaluation.findings,
            },
            "result_validation": {
                "schema_version": (
                    "chemsmart.program-result-validation-receipt.v1"
                ),
                "state": "invalid",
                "findings": evaluation.findings,
                "observations": evaluation.observations,
            },
        },
    }

    projected = project_tool_feedback(
        tool="execute_approved_program_node",
        result=tool_result,
        mode="causal-v1",
    )
    failure = projected.content["result"]["result_validation"]["observations"][
        "gaussian"
    ]["outputs"][0]["native_failure"]

    assert projected.content["status"] == "failed"
    assert failure["termination_state"] == "error_termination"
    assert failure["error_class"] == "input_syntax"
    assert failure["diagnostic_lines"] == [
        "Gaussian reported QPErr input syntax failure.",
    ]
    assert failure["artifact_refs"] == [
        {
            "artifact_id": "result.gaussian.output",
            "kind": "gaussian_output",
            "sha256": output.sha256,
            "size_bytes": output.size_bytes,
        }
    ]
