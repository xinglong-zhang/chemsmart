"""Molecule-grouper fixtures."""

import os
from unittest.mock import patch

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from tests.helpers import (
    cleanup_pymol_grouper,
    make_mock_molecules,
)


@pytest.fixture()
def make_hydrogen_molecule(hydrogen_molecule):
    """Factory for copies of ``hydrogen_molecule`` (H₂), optionally
    translated along x and/or carrying an energy.

    Grouper tests need several distinct H₂ instances per test (e.g. two
    nearly coincident ones for RMSD, or a series with different energies),
    so this returns a fresh molecule on every call.
    """

    def _make(x_offset=0.0, energy=None):
        return Molecule(
            symbols=list(hydrogen_molecule.chemical_symbols),
            positions=hydrogen_molecule.positions
            + np.array([x_offset, 0.0, 0.0]),
            energy=energy,
        )

    return _make


@pytest.fixture()
def ch_molecule(hydrogen_molecule):
    """CH on ``hydrogen_molecule``'s geometry: same atom count as H₂ but
    different elements, for element-mismatch tests."""
    return Molecule(
        symbols=["C", "H"], positions=hydrogen_molecule.positions.copy()
    )


@pytest.fixture()
def grouper_job(tmp_path, monkeypatch):
    """An RMSD ``GrouperJob`` labelled ``mygrp`` over three mock molecules,
    run from inside ``tmp_path``."""
    from chemsmart.jobs.grouper.job import GrouperJob

    monkeypatch.chdir(tmp_path)
    return GrouperJob(
        molecules=make_mock_molecules(3),
        grouping_strategy="rmsd",
        label="mygrp",
        threshold=0.5,
    )


@pytest.fixture()
def make_irmsd_grouper(make_hydrogen_molecule):
    """Factory for an ``IRMSDGrouper`` that doesn't need the irmsd binary.

    ``find_irmsd_command()`` is called at ``__init__`` time and raises
    ``RuntimeError`` if no real irmsd binary is found. Mocking it lets
    tests exercise the rest of IRMSDGrouper's logic (XYZ writing, output
    parsing, subprocess error handling, and ``_record_results``) without
    the external tool installed. Defaults to two nearly coincident H₂.
    """
    from chemsmart.jobs.grouper.rmsd import IRMSDGrouper

    def _make(molecules=None, **kwargs):
        if molecules is None:
            molecules = [
                make_hydrogen_molecule(),
                make_hydrogen_molecule(0.001),
            ]
        with patch(
            "chemsmart.jobs.grouper.rmsd.find_irmsd_command",
            return_value="/fake/bin/irmsd",
        ):
            return IRMSDGrouper(molecules, **kwargs)

    return _make


@pytest.fixture()
def make_pymol_rmsd_grouper(make_hydrogen_molecule):
    """Factory for real ``PymolRMSDGrouper`` instances, cleaned up at
    teardown (see ``cleanup_pymol_grouper``) even if the test fails.
    Defaults to two nearly coincident H₂."""
    from chemsmart.jobs.grouper.rmsd import PymolRMSDGrouper

    created = []

    def _make(molecules=None, **kwargs):
        if molecules is None:
            molecules = [
                make_hydrogen_molecule(),
                make_hydrogen_molecule(0.001),
            ]
        grouper = PymolRMSDGrouper(molecules, **kwargs)
        created.append(grouper)
        return grouper

    yield _make
    for grouper in created:
        cleanup_pymol_grouper(grouper)


@pytest.fixture()
def shared_pymol_grouper(make_pymol_rmsd_grouper):
    """One default ``PymolRMSDGrouper`` for tests that only inspect it.

    Constructing PyMOL's embedded, process-wide session is expensive and
    (empirically) prone to hangs if done many times back-to-back, so tests
    that don't need distinct construction arguments share this one.
    """
    return make_pymol_rmsd_grouper()


@pytest.fixture()
def results_recorder_args():
    """Fresh ``(header_info, sheets_data, matrix_data)`` arguments for
    ``ResultsRecorder.record_results``, with a 2x2 matrix containing
    ``inf`` off-diagonals."""
    import pandas as pd

    header_info = [
        ("Title", None),
        ("", "Some Title Line"),
        ("Key", ""),
        ("Method", "RMSD"),
    ]
    sheets_data = {"Groups": pd.DataFrame({"Group": [1], "Members": ["1"]})}
    matrix = np.array([[0.0, np.inf], [np.inf, 0.0]])
    matrix_data = ("Matrix", matrix, ["1", "2"])
    return header_info, sheets_data, matrix_data


@pytest.fixture()
def io_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "IOTests")


@pytest.fixture()
def excel_file(io_test_directory):
    return os.path.join(io_test_directory, "test.xlsx")


@pytest.fixture()
def constrained_pbc_db_file(io_test_directory):
    """Fixture of a .db file containing constrained PBC database
    from heterogeneous catalysis."""
    return os.path.join(
        io_test_directory, "heterogenous_pbc_constraints_5images.db"
    )


@pytest.fixture()
def clear_irmsd_env(monkeypatch):
    """Unset the environment variables ``find_irmsd_command`` consults
    (``IRMSD_PATH``, ``IRMSD_CONDA_ENV``, ``CONDA_PREFIX``,
    ``CONDA_PREFIX_1``), so the host's own conda/irmsd setup can't affect
    the test. Opt in with ``@pytest.mark.usefixtures("clear_irmsd_env")``."""
    for key in (
        "IRMSD_PATH",
        "IRMSD_CONDA_ENV",
        "CONDA_PREFIX",
        "CONDA_PREFIX_1",
    ):
        monkeypatch.delenv(key, raising=False)
