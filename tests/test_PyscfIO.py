"""Boundary tests for the structured PySCF output reader."""

import h5py
import numpy as np
import pytest
from ase import units

from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.writer import (
    RESULTS_SCHEMA_VERSION,
    write_pyscf_h5,
)


@pytest.fixture()
def pyscf_output(tmp_path):
    logfile = tmp_path / "water_hess.out"
    logfile.write_text("PySCF version 2.14.0\nPySCF path /opt/pyscf\n")

    spec = {
        "label": "water_hess",
        "program": "pyscf",
        "jobtype": "hess",
        "engine": "cpu",
        "stages": ["scf", "hess"],
        "symbols": ["O", "H", "H"],
        "unit": "Angstrom",
        "method": "b3lyp",
        "xc": "b3lyp",
        "basis": "def2-svp",
        "charge": 0,
        "spin": 0,
        "multiplicity": 1,
        "density_fit": True,
        "defgrid": "defgrid2",
        "solvent_model": None,
        "solvent_id": None,
        "num_basis_functions": 24,
        "num_shells": 12,
    }
    provenance = {
        "pyscf_version": "2.14.0",
        "gpu4pyscf_version": None,
        "engine": "cpu",
        "settings_digest": "b" * 64,
        "project_yaml_digest": "c" * 64,
        "wall_seconds": 2.0,
        "core_seconds": 4.0,
    }
    status = {
        "normal_termination": True,
        "stages": {
            "scf": {"converged": True, "iterations": 7},
            "hess": {"converged": True},
        },
        "failure": None,
    }
    positions = np.array(
        [
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ],
        dtype=np.float64,
    )
    results = {
        "energies": np.array([-76.3581519800948], dtype=np.float64),
        "positions": positions,
        "atomic_numbers": np.array([8, 1, 1], dtype=np.int64),
        "mo_energy": np.array([-0.50, -0.29, 0.047], dtype=np.float64),
        "mo_occ": np.array([2.0, 2.0, 0.0], dtype=np.float64),
        "vibrational_frequencies": np.array(
            [1597.96, 3905.71, 4011.86], dtype=np.float64
        ),
        "normal_modes": np.zeros((3, 3, 3), dtype=np.float64),
        "mulliken_charges": np.array([-0.29, 0.145, 0.145]),
        "dipole_moment": np.array([0.0, 0.0, 1.8]),
        "point_group": "C2v",
    }
    write_pyscf_h5(
        tmp_path / "water_hess.h5",
        spec=spec,
        provenance=provenance,
        status=status,
        results=results,
    )
    return PySCFOutput(logfile), positions


def test_v2_artifact_uses_navigable_top_level_groups(pyscf_output):
    output, _ = pyscf_output

    with h5py.File(output.resultsfile, "r") as handle:
        assert handle.attrs["schema_version"] == RESULTS_SCHEMA_VERSION
        assert set(handle) == {"provenance", "results", "spec", "status"}
        assert all(isinstance(handle[name], h5py.Group) for name in handle)
        assert handle["spec/engine"][()].decode() == "cpu"
        assert handle["status/stages/scf/iterations"][()] == 7


def test_output_boundary_preserves_documented_scientific_units(pyscf_output):
    output, positions_angstrom = pyscf_output

    assert output.energies == pytest.approx([-76.3581519800948])
    np.testing.assert_array_equal(output.positions, positions_angstrom)
    assert output.vibrational_frequencies == pytest.approx(
        [1597.96, 3905.71, 4011.86]
    )
    # Orbital-energy properties share the Gaussian/ORCA convention of eV,
    # whereas total energies remain Hartree.
    assert output.alpha_occ_eigenvalues == pytest.approx(
        [-0.50 * units.Hartree, -0.29 * units.Hartree]
    )
    assert output.alpha_virtual_eigenvalues == pytest.approx(
        [0.047 * units.Hartree]
    )


def test_output_exposes_status_and_receipt_identity(pyscf_output):
    output, _ = pyscf_output

    assert output.normal_termination is True
    assert output.failure is None
    assert output.engine == "cpu"
    assert output.version == "2.14.0"
    assert output.settings_digest == "b" * 64
    assert output.project_yaml_digest == "c" * 64
    assert output.total_elapsed_walltime == pytest.approx(2.0)
    assert output.total_core_hours == pytest.approx(4.0 / 3600.0)
    assert output.point_group == "C2v"
    assert output.rotational_symmetry_number == 2


def test_molecule_round_trip_preserves_float64_geometry(pyscf_output):
    output, positions = pyscf_output

    molecule = output.get_molecule()

    assert molecule.chemical_symbols == ["O", "H", "H"]
    assert molecule.charge == 0
    assert molecule.multiplicity == 1
    assert molecule.energy == pytest.approx(-76.3581519800948)
    assert molecule.positions.dtype == np.float64
    np.testing.assert_array_equal(molecule.positions, positions)


def test_failed_status_never_reads_as_complete(tmp_path):
    logfile = tmp_path / "failed.out"
    logfile.write_text("PySCF version 2.14.0\n")
    write_pyscf_h5(
        tmp_path / "failed.h5",
        spec={"jobtype": "sp", "stages": ["scf"]},
        provenance={"engine": "cpu"},
        status={
            "normal_termination": False,
            "stages": {"scf": {"converged": False}},
            "failure": {
                "type": "RuntimeError",
                "message": "SCF did not converge",
                "stage": "scf",
            },
        },
        results={},
    )

    output = PySCFOutput(logfile)
    assert output.normal_termination is False
    assert output.failure["stage"] == "scf"
