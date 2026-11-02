"""PySCF data-file, molecule and results fixtures."""

import os
from pathlib import Path

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.writer import (
    write_pyscf_h5,
)


# master PySCF test directory
@pytest.fixture()
def pyscf_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "PySCFTests")


# pyscf yaml files
@pytest.fixture()
def pyscf_yaml_settings_directory(pyscf_test_directory):
    return os.path.join(pyscf_test_directory, "project_yaml")


@pytest.fixture()
def pyscf_yaml_settings_test(pyscf_yaml_settings_directory):
    return os.path.join(pyscf_yaml_settings_directory, "test.yaml")


@pytest.fixture()
def pyscf_test_project_yaml(pyscf_yaml_settings_test):
    return Path(pyscf_yaml_settings_test)


@pytest.fixture()
def pyscf_inputs_directory(pyscf_test_directory):
    pyscf_inputs_directory = os.path.join(pyscf_test_directory, "inputs")
    return os.path.abspath(pyscf_inputs_directory)


@pytest.fixture()
def pyscf_outputs_directory(pyscf_test_directory):
    pyscf_outputs_directory = os.path.join(pyscf_test_directory, "outputs")
    return os.path.abspath(pyscf_outputs_directory)


@pytest.fixture()
def pyscf_water_sp_outfolder(pyscf_outputs_directory):
    return os.path.join(pyscf_outputs_directory, "water_sp_cpcm_water")


@pytest.fixture()
def pyscf_water_opt_outfolder(pyscf_outputs_directory):
    return os.path.join(pyscf_outputs_directory, "water_opt_gas_phase")


@pytest.fixture()
def pyscf_water_hess_outfolder(pyscf_outputs_directory):
    return os.path.join(pyscf_outputs_directory, "water_hess_gas_phase")


@pytest.fixture()
def pyscf_water_opt_distorted_outfolder(pyscf_outputs_directory):
    return os.path.join(pyscf_outputs_directory, "water_opt_distorted")


@pytest.fixture()
def pyscf_he_sp_outfolder(pyscf_outputs_directory):
    return os.path.join(pyscf_outputs_directory, "he_sp_gas_phase")


@pytest.fixture()
def pyscf_water_molecule(pyscf_inputs_directory):
    filepath = os.path.join(pyscf_inputs_directory, "water.xyz")
    molecule = Molecule.from_filepath(filepath)
    molecule.charge = 0
    molecule.multiplicity = 1
    return molecule


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
