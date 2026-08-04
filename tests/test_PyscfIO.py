"""Boundary tests for the structured PySCF output reader."""

import h5py
import numpy as np
import pytest
from ase import units

from chemsmart.database.assemble import PySCFAssembler
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import (
    PYSCF_SOURCE_ARTIFACT_INFO_KEY,
    PySCFOutput,
    pyscf_source_artifact_binding,
)
from chemsmart.jobs.pyscf.environment import sha256_file
from chemsmart.jobs.pyscf.writer import (
    RESULTS_SCHEMA_VERSION,
    PySCFScriptWriter,
    applied_pyscf_spec,
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
        "properties": {
            "mulliken_charges": {"status": "ok"},
            "forces": {
                "status": "unavailable",
                "failure": {
                    "type": "NotImplementedError",
                    "message": "gradient unavailable",
                },
            },
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
    assert output.beta_occ_eigenvalues == pytest.approx(
        [-0.50 * units.Hartree, -0.29 * units.Hartree]
    )
    assert output.alpha_virtual_eigenvalues == pytest.approx(
        [0.047 * units.Hartree]
    )
    assert output.beta_virtual_eigenvalues == pytest.approx(
        [0.047 * units.Hartree]
    )


def test_restricted_single_occupation_is_alpha_only(tmp_path):
    logfile = tmp_path / "hydrogen_sp.out"
    logfile.write_text("PySCF version 2.14.0\n", encoding="utf-8")
    write_pyscf_h5(
        tmp_path / "hydrogen_sp.h5",
        spec={},
        provenance={},
        status={},
        results={
            "mo_energy": np.asarray([-0.5, 0.2], dtype=np.float64),
            "mo_occ": np.asarray([1.0, 0.0], dtype=np.float64),
        },
    )

    output = PySCFOutput(logfile)

    assert output.alpha_occ_eigenvalues == pytest.approx(
        [-0.5 * units.Hartree]
    )
    assert output.beta_occ_eigenvalues == []
    assert output.alpha_virtual_eigenvalues == pytest.approx(
        [0.2 * units.Hartree]
    )
    assert output.beta_virtual_eigenvalues == pytest.approx(
        [-0.5 * units.Hartree, 0.2 * units.Hartree]
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
    assert output.property_status["mulliken_charges"]["status"] == "ok"
    assert output.property_failures == {
        "forces": {
            "type": "NotImplementedError",
            "message": "gradient unavailable",
        }
    }


def test_molecule_round_trip_preserves_float64_geometry(pyscf_output):
    output, positions = pyscf_output

    molecule = output.get_molecule()

    assert molecule.chemical_symbols == ["O", "H", "H"]
    assert molecule.charge == 0
    assert molecule.multiplicity == 1
    assert molecule.energy == pytest.approx(-76.3581519800948)
    assert molecule.positions.dtype == np.float64
    np.testing.assert_array_equal(molecule.positions, positions)


def test_direct_hdf5_and_log_entrypoints_share_authoritative_bytes(
    pyscf_output,
):
    from_log, positions = pyscf_output
    from_hdf5 = PySCFOutput(from_log.resultsfile)

    assert from_hdf5.filename == from_log.resultsfile
    assert from_hdf5.logfile == from_log.logfile
    assert from_hdf5.result_sha256 == from_log.result_sha256
    assert from_hdf5.result_sha256 == sha256_file(from_log.resultsfile)
    np.testing.assert_array_equal(from_hdf5.positions, positions)
    assert from_hdf5.spec == from_log.spec


def test_molecule_from_hdf5_carries_exact_authoritative_binding(pyscf_output):
    output, positions = pyscf_output

    molecule = Molecule.from_filepath(output.resultsfile)
    binding = molecule.info[PYSCF_SOURCE_ARTIFACT_INFO_KEY]

    np.testing.assert_array_equal(molecule.positions, positions)
    assert binding == {
        "kind": "pyscf_hdf5",
        "path": output.resultsfile,
        "sha256": sha256_file(output.resultsfile),
    }
    assert pyscf_source_artifact_binding(molecule) == binding


def test_hessian_config_binds_exact_upstream_hdf5_hash(pyscf_output):
    output, _ = pyscf_output
    molecule = Molecule.from_filepath(output.resultsfile)
    settings = type(
        "Settings",
        (),
        {
            "jobtype": "hess",
            "title": "HDF5-bound water Hessian",
            "xc": "b3lyp",
            "ab_initio": None,
            "basis": "def2-svp",
            "charge": 0,
            "multiplicity": 1,
            "defgrid": "defgrid2",
            "solvent_model": None,
            "solvent_id": None,
            "dispersion": None,
            "density_fit": False,
            "aux_basis": None,
            "scf_tol": 1.0e-9,
            "scf_maxiter": 50,
            "opt_solver": "geometric",
            "opt_maxsteps": 100,
            "engine": "cpu",
            "method_name": "b3lyp",
            "project_yaml_digest": "project-digest",
            "validate": lambda self: self,
        },
    )()
    runner = type(
        "Runner",
        (),
        {
            "num_cores": 2,
            "mem_gb": 1.0,
            "_run_id": "run-hess",
            "_run_nonce": "nonce-hess",
        },
    )()
    job = type(
        "Job",
        (),
        {
            "molecule": molecule,
            "settings": settings,
            "jobrunner": runner,
            "label": "water_hess",
            "TYPE": "pyscf_hess",
            "stages": ("scf", "hess"),
        },
    )()

    config = PySCFScriptWriter(job).build_config()
    spec = applied_pyscf_spec(config)

    assert config["input_artifact_kind"] == "pyscf_hdf5"
    assert config["input_artifact_sha256"] == sha256_file(output.resultsfile)
    assert spec["input_artifact_sha256"] == config[
        "input_artifact_sha256"
    ]


def test_database_provenance_hashes_hdf5_not_human_log(pyscf_output):
    output, _ = pyscf_output
    assembler = PySCFAssembler(filename=output.logfile)

    provenance = assembler.build_provenance()

    assert provenance["requested_source"] == output.logfile
    assert provenance["source"] == output.resultsfile
    assert provenance["source_file_hash"] == sha256_file(output.resultsfile)


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
