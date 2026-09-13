"""The program-neutral sensors are fed through the program's own reader.

``stationary_point_order_finding``, the spin observation, the basin walk
and the imaginary-mode inputs were program-neutral where they were
computed and per-program where they were fed: the ORCA branch wrote all
of them, xTB and Gaussian the count, PySCF nothing -- so a PySCF Hessian
on a saddle ended ``validated`` with no finding and no anomaly (the
archived pyscf-clean water run, and the planar-ammonia fixture below,
before 2026-09-13).  One step now reads the same facts from whatever the
reader opened.  These tests drive the verification on real PySCF bytes.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.terminal_states import (
    STATIONARY_POINT_ORDER_FINDING,
    _artifact_scan_facts,
)
from chemsmart.agent.tool_runtime import (
    HESS_STATIONARITY_GRADIENT_EH_PER_BOHR,
    CommandCompiledToolHostV1,
    _neutral_sensor_facts,
    _output_artifact_kind,
)

FIXTURES = (
    Path(__file__).resolve().parents[1] / "data" / "PySCFTests" / "outputs"
)


def _artifacts(case: str, label: str) -> tuple[TrustedArtifactRefV1, ...]:
    """Every file the run left, typed the way the executor types them."""

    found = []
    for path in sorted((FIXTURES / case).glob(f"{label}*")):
        if path.suffix == ".json" and "reference" in path.name:
            continue
        kind = (
            "pyscf_hdf5"
            if path.suffix == ".h5"
            else _output_artifact_kind("pyscf", path)
        )
        found.append(
            TrustedArtifactRefV1(
                artifact_id=f"result.{path.name}",
                kind=kind,
                sha256=file_sha256(path),
                size_bytes=path.stat().st_size,
                path=str(path),
                cli_value=str(path),
            )
        )
    return tuple(found)


def _evaluate(
    case: str, label: str, *, jobtype: str, charge=0, multiplicity=1
):
    return CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="pyscf",
        jobtype=jobtype,
        charge=charge,
        multiplicity=multiplicity,
        output_artifacts=_artifacts(case, label),
        exit_status=0,
    )


def _signals(evaluation) -> dict[str, dict]:
    return {item["signal_id"]: item for item in evaluation.anomalies}


@pytest.mark.capability("signal:stationary_point.unexpected_order")
@pytest.mark.capability("gate:stationary_point_order")
def test_a_pyscf_hessian_on_a_saddle_is_typed_and_its_mode_named():
    evaluation = _evaluate(
        "nh3_planar_hess", "nh3_planar_hess_gas_phase", jobtype="hess"
    )
    assert STATIONARY_POINT_ORDER_FINDING in evaluation.findings
    block = evaluation.observations["pyscf"]
    assert block["consequential_imaginary_mode_count"] == 1
    assert block["imaginary_frequencies_cm1"][0] < -800.0
    assert block["charge"] == 0 and block["multiplicity"] == 1
    anomaly = _signals(evaluation)["stationary_point.unexpected_order"]
    assert anomaly["expected_imaginary_modes"] == 0
    assert anomaly["observed_imaginary_modes"] == 1
    assert anomaly["lowest_imaginary_cm1"] < -800.0
    assert "heavy_atom_share" in anomaly


@pytest.mark.capability("gate:stationary_point_order")
def test_a_pyscf_hessian_on_a_minimum_is_not_typed():
    evaluation = _evaluate(
        "water_hess", "water_hess_gas_phase", jobtype="hess"
    )
    assert STATIONARY_POINT_ORDER_FINDING not in evaluation.findings
    assert (
        evaluation.observations["pyscf"]["consequential_imaginary_mode_count"]
        == 0
    )
    signals = _signals(evaluation)
    assert "stationary_point.unexpected_order" not in signals
    assert "stationary_point.gradient_above_optimizer_criterion" not in signals


@pytest.mark.capability(
    "signal:stationary_point.gradient_above_optimizer_criterion"
)
@pytest.mark.capability("policy:hess_stationarity_gradient")
def test_zero_imaginary_modes_at_a_non_stationary_geometry_is_an_observation():
    evaluation = _evaluate(
        "water_stretched_hess",
        "water_stretched_hess_gas_phase",
        jobtype="hess",
    )
    assert (
        STATIONARY_POINT_ORDER_FINDING not in evaluation.findings
    ), "three real frequencies: the order rule has nothing to say"
    anomaly = _signals(evaluation)[
        "stationary_point.gradient_above_optimizer_criterion"
    ]
    assert anomaly["max_abs_gradient_eh_per_bohr"] > 0.018
    assert (
        anomaly["optimizer_criterion_eh_per_bohr"]
        == HESS_STATIONARITY_GRADIENT_EH_PER_BOHR
        == 4.5e-4
    )
    assert anomaly["policy_id"] == "hess_stationarity_gradient"
    assert (
        evaluation.observations["pyscf"]["max_abs_gradient_eh_per_bohr"]
        > 0.018
    )


@pytest.mark.capability("signal:spin.s2_deviation_ge_0.2")
def test_a_pyscf_radical_carries_its_spin_observation():
    evaluation = _evaluate(
        "hydroxyl_sp", "hydroxyl_sp_gas_phase", jobtype="sp", multiplicity=2
    )
    block = evaluation.observations["pyscf"]
    assert abs(block["spin_square_observed"] - 0.7518) < 1e-3
    assert block["spin_square_expected"] == 0.75
    assert abs(block["spin_square_deviation"]) < 0.01
    assert "spin.s2_deviation_ge_0.2" not in _signals(evaluation)


@pytest.mark.capability("signal:geometry.connectivity_changed")
def test_the_basin_sensor_reaches_a_pyscf_optimisation():
    """Measured against the geometry the node was handed, through the
    same helper the ORCA branch used to be the only caller of."""

    xyz = FIXTURES / "inputs" / "water_distorted.xyz"
    supplied = TrustedArtifactRefV1(
        artifact_id="geometry.water.distorted",
        kind="geometry_xyz",
        sha256=file_sha256(xyz),
        size_bytes=xyz.stat().st_size,
        path=str(xyz),
        cli_value=str(xyz),
    )
    block, inputs = _neutral_sensor_facts(
        program="pyscf",
        jobtype="opt",
        multiplicity=1,
        output_artifacts=_artifacts("water_opt", "water_opt_gas_phase"),
        expected_input_artifact=supplied,
        expected_root_artifact=None,
    )
    assert block["charge"] == 0 and block["multiplicity"] == 1
    assert inputs["basin"]["bonds_made"] == []
    assert (
        inputs["basin"]["bonds_broken"] == []
    ), "a 1.10 A / 90 deg water relaxing to 0.967 A keeps its two bonds"


def test_convergence_facts_are_read_through_the_readers():
    def record(case, name):
        path = FIXTURES / case / name
        return {
            "program": "pyscf",
            "output_artifacts": [
                {
                    "kind": "pyscf_hdf5",
                    "path": str(path),
                    "sha256": file_sha256(path),
                }
            ],
        }

    converged, _reached, _planned, digests = _artifact_scan_facts(
        record("water_opt", "water_opt_gas_phase.h5")
    )
    assert converged is True and len(digests) == 1
    converged, *_ = _artifact_scan_facts(
        record("water_opt_maxsteps1", "water_opt_maxsteps1_gas_phase.h5")
    )
    assert converged is False
    converged, *_ = _artifact_scan_facts(
        record("water_sp", "water_sp_gas_phase.h5")
    )
    assert converged is None, "a single point has no optimisation to converge"


# ----------------------------------------------------------------------
# the response stage's numbers are sensor facts: no threshold, no signal
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:opt:excited_state_followed_root")
def test_an_excited_root_optimisation_reports_its_gaps_as_facts():
    """Which root was followed, how far it ended from the ground state and
    from its neighbour, and how many roots the filter dropped -- read
    through the reader into the observation block, with no policy behind
    them.  A root is an index, never a state identity; the session reads
    the gap before it calls the delivered geometry a state's minimum."""

    evaluation = _evaluate(
        "formaldehyde_s1_opt_planar",
        "formaldehyde_s1_opt_planar_gas_phase",
        jobtype="opt",
    )
    block = evaluation.observations["pyscf"]
    assert block["excited_state_followed_root"] == 1
    assert abs(block["excited_state_root_gap_to_ground_ev"] - 3.406) < 5e-3
    assert block["excited_state_root_gap_to_neighbour_ev"] > 4.0
    assert block["excited_state_roots_filtered"] == 0
    assert block["excited_state_followed_root_converged"] is True
    assert block["excited_state_final_gradient_max"] < 4.5e-4
    assert not [
        item for item in evaluation.anomalies if "excited" in item["signal_id"]
    ], "the gap is a fact, not a signal: no policy exists for it"

    degenerate = _evaluate(
        "water_s1_opt_degenerate",
        "water_s1_opt_degenerate_gas_phase",
        jobtype="opt",
    ).observations["pyscf"]
    assert degenerate["excited_state_root_gap_to_neighbour_ev"] < 1e-4


@pytest.mark.capability("selector:pyscf:td:excitation_energies")
def test_a_spectrum_reports_its_lowest_root_and_the_filter_count():
    block = _evaluate(
        "water_td_singlet", "water_td_singlet_gas_phase", jobtype="td"
    ).observations["pyscf"]
    assert abs(block["excited_state_lowest_root_ev"] - 7.552) < 5e-3
    assert block["excited_state_roots_requested"] == 3
    assert block["excited_state_roots_obtained"] == 3
    assert block["excited_state_roots_filtered"] == 0
    assert block["excited_state_unconverged_roots"] == []
    assert "excited_state_followed_root" not in block
    ground = _evaluate(
        "water_opt", "water_opt_gas_phase", jobtype="opt"
    ).observations["pyscf"]
    assert not [key for key in ground if key.startswith("excited_state")]
