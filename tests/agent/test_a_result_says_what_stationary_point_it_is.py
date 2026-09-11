"""A failed search keeps its failure and can still say what it found.

Twenty-four archived optimisations ended on a structural saddle and none
was delivered as the stationary point it was. The numbers were always
readable -- a wrong-stationary-point run terminates normally -- so what
was missing was not access but a way for the delivered number to say
which structure it belongs to. The session states the order; the host
checks it against the frequencies the program printed, under the same
20 cm^-1 convention the validator uses, and mints a receipt beside the
failure rather than in place of it.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.agent.execution import (
    build_stationary_point_characterisation,
)

_OUT = Path(__file__).resolve().parents[1] / "data" / "ORCATests" / "outputs"


def _artifact(name: str, kind: str = "orca_output") -> TrustedArtifactRefV1:
    path = _OUT / name
    return TrustedArtifactRefV1(
        artifact_id=f"result.{name}",
        kind=kind,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


@pytest.mark.capability("tool:characterise_stationary_point")
def test_a_saddle_may_say_it_is_a_saddle():
    receipt = build_stationary_point_characterisation(
        result_artifact=_artifact("sn2_ts.out"),
        program="orca",
        order_claimed=1,
        node_id="opt-a",
    )
    assert receipt.observed_imaginary_modes == 1
    assert receipt.order_claimed == 1
    assert receipt.lowest_imaginary_cm_1 == pytest.approx(-407.58)
    assert receipt.node_id == "opt-a"
    assert receipt.receipt_sha256


@pytest.mark.capability("tool:characterise_stationary_point")
def test_the_printed_modes_refuse_an_order_they_do_not_support():
    with pytest.raises(ContractError, match="1 imaginary mode"):
        build_stationary_point_characterisation(
            result_artifact=_artifact("sn2_ts.out"),
            program="orca",
            order_claimed=0,
        )


@pytest.mark.capability("tool:characterise_stationary_point")
def test_a_true_minimum_characterises_as_one():
    receipt = build_stationary_point_characterisation(
        result_artifact=_artifact("water_opt.out"),
        program="orca",
        order_claimed=0,
    )
    assert receipt.observed_imaginary_modes == 0
    assert receipt.lowest_imaginary_cm_1 is None


@pytest.mark.capability("tool:characterise_stationary_point")
def test_a_result_with_no_frequencies_characterises_nothing():
    with pytest.raises(ContractError, match="printed no frequencies"):
        build_stationary_point_characterisation(
            result_artifact=_artifact("fe3_sextet.out"),
            program="orca",
            order_claimed=0,
        )


@pytest.mark.capability("tool:characterise_stationary_point")
def test_the_receipt_is_digest_bound_and_names_its_result():
    from dataclasses import replace

    receipt = build_stationary_point_characterisation(
        result_artifact=_artifact("sn2_ts.out"),
        program="orca",
        order_claimed=1,
    )
    assert receipt.result_artifact_sha256 == _artifact("sn2_ts.out").sha256
    with pytest.raises(ContractError, match="digest mismatch"):
        replace(receipt, order_claimed=2)


@pytest.mark.capability("tool:characterise_stationary_point")
def test_it_moves_no_verdict(tmp_path):
    """The node's own validation is byte-identical either way: the
    characterisation stands beside the failure, never in place of it."""

    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    artifact = _artifact("sn2_ts.out")
    before = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="opt",
        charge=-1,
        multiplicity=1,
        output_artifacts=(artifact,),
        exit_status=0,
    )
    build_stationary_point_characterisation(
        result_artifact=artifact, program="orca", order_claimed=1
    )
    after = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="opt",
        charge=-1,
        multiplicity=1,
        output_artifacts=(artifact,),
        exit_status=0,
    )
    assert before.findings == after.findings
    assert "result.stationary_point_order" in before.findings
