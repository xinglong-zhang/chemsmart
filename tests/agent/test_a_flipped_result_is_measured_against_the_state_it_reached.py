"""<S**2> of a flipped result is measured against the Ms it converged to.

A flip of named centres (``site_spin_flip``, ORCA's FlipSpin with FinalMs,
or a native BrokenSym) starts from the high-spin determinant the coordinate
line names and converges to the Ms it asks for; ORCA prints which ("converge
to the broken symmetry state with Ms= 0.0"). Such a determinant is an
eigenfunction of S_z at that Ms and a mixture of total spins, and the pure
state it stands for is S = |Ms|: for two S = 1 centres coupled
antiferromagnetically, the singlet. The host measured its <S**2> against
the coordinate line's S(S+1) instead -- 6.0 for the quintet two Ni(II) start
from -- and raised ``spin.s2_deviation_ge_0.2`` at -4.0 on every flipped
result, an anomaly about a state the run never targeted (R11 E2's typed
flip; ino2's Ni(II)2, CUHK 2157086 and 2157103). It also did not count the
flip as the broken-symmetry request it is, so a flip that collapsed to a
spin-symmetric solution could never raise
``spin.broken_symmetry_request_unbroken``.

Measured against S = |Ms| a flipped Ni(II)2 reads as a GuessMix singlet
does: ``broken``, with the contamination a two-site broken-symmetry
determinant carries by construction (<S**2> of about S_A + S_B = 2, here
1.998). The fixtures are the oracle's real ORCA 6.1.1 outputs; the controls
are a GuessMix singlet and an ordinary sextet, whose readings do not move.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.analysis.result_readers import reader_for

_DATA = Path(__file__).resolve().parents[1] / "data" / "ORCATests"
_BS = _DATA / "broken_symmetry"

pytestmark = pytest.mark.capability("signal:spin.s2_deviation_ge_0.2")


def _evaluate(path: Path, *, multiplicity: int, charge: int):
    artifact = TrustedArtifactRefV1(
        artifact_id="result.spin",
        kind="orca_output",
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )
    return CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="sp",
        charge=charge,
        multiplicity=multiplicity,
        output_artifacts=(artifact,),
        exit_status=0,
    )


def _signals(evaluation, signal_id):
    return [
        item for item in evaluation.anomalies if item["signal_id"] == signal_id
    ]


@pytest.mark.parametrize(
    "name", ["ni2_fs1_gas_phase.out", "ni2_bs22_gas_phase.out"]
)
def test_a_flipped_result_is_measured_against_the_ms_it_converged_to(name):
    evaluation = _evaluate(_BS / name, multiplicity=5, charge=2)
    observation = evaluation.observations["orca"]
    # The pure state an Ms = 0 determinant stands for is S = 0.
    assert observation["spin_square_expected"] == 0.0
    assert observation["spin_square_deviation"] == pytest.approx(
        1.998, abs=1e-3
    )
    (anomaly,) = _signals(evaluation, "spin.s2_deviation_ge_0.2")
    assert anomaly["spin_square_deviation"] == pytest.approx(1.998, abs=1e-3)
    # The record says what the deviation was measured against.
    assert anomaly["final_ms"] == 0.0
    assert anomaly["bound_multiplicity"] == 5
    # A flip is the broken-symmetry request it is, and it broke.
    output = reader_for("orca").open_output(str(_BS / name))
    level = reader_for("orca").level_for_output(output)
    assert level["broken_symmetry"] is True
    assert level["final_ms"] == 0.0
    assert not _signals(evaluation, "spin.broken_symmetry_request_unbroken")


def test_the_typed_flip_is_named_on_the_level_it_ran_at():
    """What ran, in the host's own atom numbering: the second Ni."""

    reader = reader_for("orca")
    output = reader.open_output(str(_BS / "ni2_fs1_gas_phase.out"))
    level = reader.level_for_output(output)
    assert level["site_spin_flip"] == {"atoms": [2], "final_ms": 0.0}


@pytest.mark.parametrize(
    "path, multiplicity, charge, expected, deviation",
    [
        # A GuessMix singlet: its target was always S = 0.
        (_BS / "o_pbenzyne_bs_gas_phase.out", 1, 0, 0.0, 0.970279),
        # An ordinary sextet: S = 5/2, as the coordinate line says.
        (_DATA / "outputs" / "fe3_sextet.out", 6, 3, 8.75, 0.009007),
    ],
)
def test_an_unflipped_result_keeps_its_coordinate_line_target(
    path, multiplicity, charge, expected, deviation
):
    evaluation = _evaluate(path, multiplicity=multiplicity, charge=charge)
    observation = evaluation.observations["orca"]
    assert observation["spin_square_expected"] == pytest.approx(expected)
    assert observation["spin_square_deviation"] == pytest.approx(
        deviation, abs=1e-6
    )
    assert "final_ms" not in reader_for("orca").level_for_output(
        reader_for("orca").open_output(str(path))
    )
