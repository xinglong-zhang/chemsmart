"""xtb data-file and input fixtures."""

import copy
import os

import pytest

from tests.helpers import bare_xtb_input


# master xTB test directory
@pytest.fixture()
def xtb_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "XTBTests")


@pytest.fixture()
def xtb_inputs_directory(xtb_test_directory):
    xtb_inputs_directory = os.path.join(xtb_test_directory, "inputs")
    return os.path.abspath(xtb_inputs_directory)


@pytest.fixture()
def xtb_default_inputfile(xtb_inputs_directory):
    return os.path.join(xtb_inputs_directory, "default.inp")


@pytest.fixture()
def xtb_sp_alpb_inputfile(xtb_inputs_directory):
    return os.path.join(xtb_inputs_directory, "alpb_water.inp")


@pytest.fixture()
def xtb_outputs_directory(xtb_test_directory):
    xtb_outputs_directory = os.path.join(xtb_test_directory, "outputs")
    return os.path.abspath(xtb_outputs_directory)


@pytest.fixture()
def xtb_co2_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "co2_ohess")


@pytest.fixture()
def xtb_water_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "water_ohess")


@pytest.fixture()
def xtb_cyclopentadienyl_anion_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "cyclopentadienyl_anion_opt")


@pytest.fixture()
def xtb_p_benzyne_opt_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "p_benzyne_opt_alpb_toluene")


@pytest.fixture()
def xtb_p_benzyne_sp_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "p_benzyne_sp_alpb_toluene")


@pytest.fixture()
def xtb_acetaldehyde_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "acetaldehyde_hess")


@pytest.fixture()
def xtb_he_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "he_hess")


_MALFORMED_GROUPS = {
    "chrg": ["notanumber"],
    "spin": ["notanumber"],
    "cube": [
        "step=notanumber",
        "pthr=notanumber",
        "boff=notanumber",
        "cal=notanumber",
    ],
    "embedding": ["at=notanumber"],
    "gfn": ["dispscale=notanumber"],
    "hess": [
        "sccacc=notanumber",
        "step=notanumber",
        "scale=notanumber",
    ],
    "md": [
        "temp=notanumber",
        "time=notanumber",
        "dump=notanumber",
        "skip=notanumber",
        "step=notanumber",
        "hmass=notanumber",
        "shake=notanumber",
        "sccacc=notanumber",
    ],
    "modef": [
        "n=notanumber",
        "step=notanumber",
        "updat=notanumber",
        "local=notanumber",
        "vthr=notanumber",
        "prj=notanumber",
        "mode=notanumber",
    ],
    "opt": [
        "optlevel=notanumber",
        "microcycle=notanumber",
        "maxcycle=notanumber",
        "maxdispl=notanumber",
        "hlow=notanumber",
        "s6=notanumber",
        "kstretch=notanumber",
        "kbend=notanumber",
        "ktorsion=notanumber",
        "koutofp=notanumber",
        "kvdw=notanumber",
        "kes=notanumber",
        "rcut=notanumber",
    ],
    "path": [
        "nrun=notanumber",
        "npoint=notanumber",
        "anopt=notanumber",
        "kpush=notanumber",
        "kpull=notanumber",
        "alp=notanumber",
    ],
    "scc": [
        "maxiterations=notanumber",
        "temp=notanumber",
        "broydamp=notanumber",
    ],
    "symmetry": ["desy=notanumber", "maxat=notanumber"],
    "thermo": [
        "temp=notanumber",
        "sthr=notanumber",
        "imagthr=notanumber",
        "scale=notanumber",
    ],
    "wall": [
        "alpha=notanumber",
        "beta=notanumber",
        "temp=notanumber",
        "autoscale=notanumber",
        "axisshift=notanumber",
    ],
}


@pytest.fixture()
def malformed_xtb_input():
    """Bare XTBInput with every numeric key set to a non-numeric
    string, to force every ValueError-fallback branch at once."""
    return bare_xtb_input(copy.deepcopy(_MALFORMED_GROUPS))
