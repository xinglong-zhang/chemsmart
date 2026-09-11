"""<S²> is recorded beside the bound state, and a large deviation has standing.

The charter said the spin expectation value is recorded as an observation
beside the bound multiplicity. It was not: the readers expose the table as
a property, the host called it as a method, and the TypeError was swallowed,
so the observation had been empty on every ORCA and Gaussian result. Now
the value is read, and a deviation of 0.2 or more from S(S+1) is an anomaly
observation with the number that tripped it -- never a finding, never a
gate, because a requested broken-symmetry state is not a surprise.
"""

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_DATA = Path(__file__).resolve().parents[1] / "data"
_ORCA = _DATA / "ORCATests" / "outputs"
_GAUSSIAN = _DATA / "GaussianTests" / "outputs"


def _evaluate(
    path: Path,
    multiplicity: int,
    *,
    program: str = "orca",
    jobtype: str = "opt",
    charge: int = 3,
):
    artifact = TrustedArtifactRefV1(
        artifact_id="result.iron",
        kind=f"{program}_output",
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )
    return CommandCompiledToolHostV1._evaluate_execution_outputs(
        program=program,
        jobtype=jobtype,
        charge=charge,
        multiplicity=multiplicity,
        output_artifacts=(artifact,),
        exit_status=0,
    )


def _spin_anomalies(evaluation):
    return [
        item
        for item in evaluation.anomalies
        if item["signal_id"] == "spin.s2_deviation_ge_0.2"
    ]


@pytest.mark.capability("signal:spin.s2_deviation_ge_0.2")
def test_the_spin_expectation_value_is_recorded_beside_the_bound_state():
    observation = _evaluate(_ORCA / "fe3_sextet.out", 6).observations["orca"]
    assert observation["spin_square_expected"] == pytest.approx(8.75)
    assert observation["spin_square_observed"] == pytest.approx(8.759007)
    assert observation["spin_square_deviation"] == pytest.approx(
        0.009007, abs=1e-6
    )


@pytest.mark.capability("signal:spin.s2_deviation_ge_0.2")
def test_a_large_deviation_is_an_observation_with_its_number():
    evaluation = _evaluate(_ORCA / "fe3_doublet.out", 2)
    (anomaly,) = _spin_anomalies(evaluation)
    assert anomaly["spin_square_deviation"] == pytest.approx(0.9547)
    assert anomaly["bound_multiplicity"] == 2
    # The deviation moves no finding: it is weighed, never judged.
    assert not any("spin" in finding for finding in evaluation.findings)


@pytest.mark.capability("signal:spin.s2_deviation_ge_0.2")
def test_a_small_deviation_is_no_anomaly():
    assert _spin_anomalies(_evaluate(_ORCA / "fe3_sextet.out", 6)) == []


@pytest.mark.capability("signal:spin.s2_deviation_ge_0.2")
def test_the_sensor_reads_every_program_whose_spin_the_host_reads():
    """Gaussian records the value before and after annihilation and
    keeps it per output row; the sensor reads the wavefunction's own
    value there as it reads ORCA's, or the rule is a hole."""

    evaluation = _evaluate(
        _GAUSSIAN / "iron_neutral_triplet.log",
        3,
        program="gaussian",
        jobtype="ts",
        charge=0,
    )
    (anomaly,) = _spin_anomalies(evaluation)
    assert anomaly["spin_square_deviation"] == pytest.approx(0.4589)
    assert anomaly["bound_multiplicity"] == 3
