"""PyMOL fixtures."""

import shutil

import pytest


############ PyMOL Fixtures ##################
@pytest.fixture(scope="session")
def skip_if_no_pymol():
    if shutil.which("pymol") is None:
        pytest.skip("PyMOL not installed")


@pytest.fixture()
def loaded_1_mer(visualized_1_mer_xyz_file):
    """Load ``1-mer.xyz`` into a fresh, real PyMOL session (PyMOL required)."""
    from pymol import cmd

    cmd.reinitialize()
    cmd.load(visualized_1_mer_xyz_file, "mol1mer")
    yield "mol1mer"
    cmd.delete("all")
