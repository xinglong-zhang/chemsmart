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


def _opt_status(optimizer_converged):
    return {
        "normal_termination": optimizer_converged,
        "stages": {
            "scf": {"converged": True},
            "opt": {
                "converged": optimizer_converged,
                "optimizer_converged": optimizer_converged,
                "final_scf_converged": True,
            },
        },
        "failure": None,
    }


def test_unconverged_optimisation_is_not_labelled_optimized(tmp_path):
    logfile = tmp_path / "opt.out"
    logfile.write_text("PySCF version 2.14.0\n")
    write_pyscf_h5(
        tmp_path / "opt.h5",
        spec={
            "jobtype": "opt",
            "symbols": ["H", "H"],
            "charge": 0,
            "multiplicity": 1,
        },
        provenance={"engine": "cpu"},
        status=_opt_status(False),
        results={"positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]},
    )

    output = PySCFOutput(logfile)
    assert output.get_molecule().is_optimized_structure is False
    assert output.optimized_structure is None


def test_run_without_a_geometry_reports_why(tmp_path):
    logfile = tmp_path / "failed.out"
    logfile.write_text("PySCF version 2.14.0\n")
    write_pyscf_h5(
        tmp_path / "failed.h5",
        spec={"jobtype": "sp", "symbols": ["H", "H"]},
        provenance={"engine": "cpu"},
        status={
            "normal_termination": False,
            "stages": {},
            "failure": {
                "type": "RuntimeError",
                "message": "SMD module is not available",
                "stage": "initialization",
            },
        },
        results={},
    )

    with pytest.raises(ValueError, match="SMD module is not available"):
        PySCFOutput(logfile).get_molecule()
