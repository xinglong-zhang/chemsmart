"""Boundary tests for the structured PySCF output reader."""

import logging

import h5py
import numpy as np
import pytest
from ase import units

from chemsmart.io.pyscf.output import (
    PySCFOutput,
    read_pyscf_h5,
    rotational_symmetry_number_from_point_group,
)
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


def test_missing_results_read_as_not_terminated(tmp_path):
    logfile = tmp_path / "absent.out"
    logfile.write_text("PySCF version 2.14.0\n")
    output = PySCFOutput(logfile)
    assert output.normal_termination is False
    assert output.file_date is None
    with pytest.raises(FileNotFoundError, match="did not reach"):
        output.spec


def test_corrupt_results_read_as_not_terminated(tmp_path):
    """A file truncated by a crash must never count as a completed run."""
    logfile = tmp_path / "corrupt.out"
    logfile.write_text("PySCF version 2.14.0\n")
    (tmp_path / "corrupt.h5").write_bytes(b"not an hdf5 file")
    assert PySCFOutput(logfile).normal_termination is False


def test_unsupported_schema_version_is_rejected(tmp_path):
    path = tmp_path / "future.h5"
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = "99.0"
    with pytest.raises(ValueError, match="Unsupported PySCF results schema"):
        read_pyscf_h5(path)


def test_missing_schema_version_is_rejected(tmp_path):
    path = tmp_path / "unversioned.h5"
    with h5py.File(path, "w"):
        pass
    with pytest.raises(ValueError, match="Unsupported PySCF results schema"):
        read_pyscf_h5(path)


def test_missing_section_is_named(tmp_path):
    path = tmp_path / "partial.h5"
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = RESULTS_SCHEMA_VERSION
        handle.create_group("spec")
    with pytest.raises(ValueError, match="provenance"):
        read_pyscf_h5(path)


@pytest.mark.parametrize(
    ("point_group", "sigma"),
    [
        ("C1", 1),
        ("Cs", 1),
        ("Ci", 1),
        ("C2", 2),
        ("C2v", 2),
        ("C3h", 3),
        ("D2", 4),
        ("D3d", 6),
        ("D6h", 12),
        ("S4", 2),
        ("S6", 3),
        ("Coov", 1),
        ("Dooh", 2),
        ("T", 12),
        ("Td", 12),
        ("Th", 12),
        ("O", 24),
        ("Oh", 24),
        ("I", 60),
        ("Ih", 60),
    ],
)
def test_rotational_symmetry_numbers(point_group, sigma):
    assert rotational_symmetry_number_from_point_group(point_group) == sigma


@pytest.mark.parametrize("point_group", [None, "", "Xyz", "Q5"])
def test_unknown_point_group_warns_and_uses_one(point_group, caplog):
    with caplog.at_level(logging.WARNING):
        sigma = rotational_symmetry_number_from_point_group(point_group)
    assert sigma == 1
    assert "Rotational entropy may be overestimated" in caplog.text
