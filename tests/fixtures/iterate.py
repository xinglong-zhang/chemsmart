"""Iterate job fixtures."""

import importlib
import os
from unittest.mock import MagicMock

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule


############ Iterate Fixtures ##################
@pytest.fixture()
def iterate_test_directory(test_data_directory):
    """Returns the absolute path to tests/data/IterateTests."""
    return os.path.join(test_data_directory, "IterateTests")


@pytest.fixture()
def iterate_input_directory(iterate_test_directory):
    """Returns the absolute path to tests/data/IterateTests/input."""
    return os.path.join(iterate_test_directory, "input")


@pytest.fixture()
def iterate_expected_output_directory(iterate_test_directory):
    """Returns the absolute path to tests/data/IterateTests/expected_output."""
    return os.path.join(iterate_test_directory, "expected_output")


@pytest.fixture()
def iterate_configs_directory(iterate_test_directory):
    """Returns the absolute path to tests/data/IterateTests/configs."""
    return os.path.join(iterate_test_directory, "configs")


@pytest.fixture()
def iterate_integration_config_file(iterate_configs_directory):
    """Returns the absolute path to
    tests/data/IterateTests/configs/integration_iterate.toml."""
    return os.path.join(iterate_configs_directory, "integration_iterate.toml")


@pytest.fixture()
def iterate_timeout_config_file(iterate_configs_directory):
    """Returns the absolute path to
    tests/data/IterateTests/configs/timeout_iterate.toml."""
    return os.path.join(iterate_configs_directory, "timeout_iterate.toml")


@pytest.fixture()
def iterate_template_file(iterate_configs_directory):
    """Returns the absolute path to
    tests/data/IterateTests/configs/iterate_template.toml."""
    return os.path.join(iterate_configs_directory, "iterate_template.toml")


@pytest.fixture()
def iterate_invalid_skeleton_link_index_config_file(iterate_configs_directory):
    """Returns the absolute path to tests/data/IterateTests/configs/
    invalid_skeleton_link_index.toml."""
    return os.path.join(
        iterate_configs_directory, "invalid_skeleton_link_index.toml"
    )


@pytest.fixture()
def iterate_expected_output_file(iterate_expected_output_directory):
    """Returns the absolute path to tests/data/IterateTests/expected_output/
    integration_iterate_SLSQP_lagrange_multipliers_96_6.xyz."""
    return os.path.join(
        iterate_expected_output_directory,
        "integration_iterate_SLSQP_lagrange_multipliers_96_6.xyz",
    )


# ── Iterate molecules (preprocessor / analyzer unit tests) ──
@pytest.fixture()
def ch2_fragment_molecule():
    """CH₂ fragment: the carbon has 2 of its 4 bonds filled, so it has an
    available bonding position (preprocessors leave it unchanged)."""
    return Molecule(
        symbols=["C", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [1.09, 0.0, 0.0],
            [-0.36, 1.03, 0.0],
        ],
    )


@pytest.fixture()
def carbon_ring_molecule():
    """Three-membered carbon ring (C₃ backbone, no H): removing any single
    bond never disconnects the ring."""
    return Molecule(
        symbols=["C", "C", "C"],
        positions=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.51, 0.0, 0.0],
                [0.755, 1.308, 0.0],
            ]
        ),
    )


@pytest.fixture()
def isolated_carbon_molecule():
    """C-H pair plus a far-away carbon (index 2) that has no bonds."""
    return Molecule(
        symbols=["C", "H", "C"],
        positions=[
            [0.0, 0.0, 0.0],
            [1.09, 0.0, 0.0],
            [20.0, 20.0, 20.0],
        ],
    )


@pytest.fixture()
def frozen_coh_molecule():
    """Linear C-O-H fragment with ``frozen_atoms=[1, 0, 1]``."""
    return Molecule(
        symbols=["C", "O", "H"],
        positions=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.4, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
            ]
        ),
        frozen_atoms=[1, 0, 1],
    )


@pytest.fixture()
def frozen_hydrogen_molecule(hydrogen_molecule):
    """``hydrogen_molecule`` with only its first atom frozen."""
    return Molecule(
        symbols=list(hydrogen_molecule.chemical_symbols),
        positions=hydrogen_molecule.positions.copy(),
        frozen_atoms=[1, 0],
    )


@pytest.fixture()
def fluorine_atom_molecule():
    """A lone F atom at the origin: the smallest possible substituent."""
    return Molecule(symbols=["F"], positions=[[0.0, 0.0, 0.0]])


# ── Iterate analyzer (Lagrange optimisation) ──
@pytest.fixture()
def iterate_skeleton_array():
    """Ethylene-like skeleton as an ``(N, 4)`` ``[Z, x, y, z]`` array, the
    layout ``IterateAnalyzer._find_optimal_position`` works on."""
    return np.array(
        [
            [6.0, 0.0, 0.0, 0.0],
            [6.0, 1.5, 0.0, 0.0],
            [1.0, -0.5, 0.9, 0.0],
            [1.0, -0.5, -0.9, 0.0],
        ]
    )


@pytest.fixture()
def iterate_single_atom_skeleton_array():
    """One-atom ``[Z, x, y, z]`` skeleton: its centre of mass coincides
    with the link atom."""
    return np.array([[6.0, 0.0, 0.0, 0.0]])


@pytest.fixture()
def methanol_fluorine_iterate_analyzer(
    methanol_molecule, fluorine_atom_molecule
):
    """Cheap ``IterateAnalyzer`` placing an F atom on methanol's carbon,
    with minimal sphere/rotation sampling."""
    from chemsmart.jobs.iterate.iterate import IterateAnalyzer

    return IterateAnalyzer(
        skeleton=methanol_molecule,
        substituent=fluorine_atom_molecule,
        skeleton_link_index=1,
        substituent_link_index=1,
        sphere_direction_samples_num=4,
        axial_rotations_sample_num=2,
    )


@pytest.fixture()
def failed_minimize_result():
    """Stand-in for an unsuccessful ``scipy.optimize.minimize`` result."""
    from types import SimpleNamespace

    return SimpleNamespace(
        success=False,
        message="did not converge",
        fun=float("inf"),
        x=np.zeros(6),
    )


@pytest.fixture()
def mock_iterate_molecule_pool():
    """``IterateMoleculePool`` holding one mock skeleton and one mock
    substituent, each returning itself from ``copy()``."""
    from chemsmart.jobs.iterate.runner import IterateMoleculePool

    skeleton = MagicMock(name="skeleton")
    skeleton.copy.return_value = skeleton
    substituent = MagicMock(name="substituent")
    substituent.copy.return_value = substituent
    return IterateMoleculePool(
        skeletons=[skeleton], substituents=[substituent]
    )


@pytest.fixture()
def patch_iterate_multiprocessing(monkeypatch):
    """Route ``IterateJobRunner.run_combinations``' ``Manager().Queue()``
    and ``Process(...)`` to the given fakes for the rest of the test:
    ``patch_iterate_multiprocessing(fake_queue, fake_process)``."""

    def _patch(fake_queue, fake_process):
        fake_manager = MagicMock()
        fake_manager.Queue.return_value = fake_queue
        fake_manager.shutdown.return_value = None
        monkeypatch.setattr(
            "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
            MagicMock(return_value=fake_manager),
        )
        monkeypatch.setattr(
            "chemsmart.jobs.iterate.runner.multiprocessing.Process",
            MagicMock(return_value=fake_process),
        )

    return _patch


# ── Iterate CLI ──
@pytest.fixture()
def patch_iterate_cli_job(mocker):
    """Mock the job pipeline behind ``chemsmart run iterate`` (job, runner,
    settings, ``validate_config``) and return the mock ``IterateJob``
    whose ``run()`` returns ``run_return``:
    ``mock_job = patch_iterate_cli_job(run_return="out.xyz")``.

    The ``chemsmart.cli.iterate`` package rebinds ``iterate`` to the Click
    group, so the submodule is looked up via ``importlib`` and patched with
    ``patch.object`` rather than by dotted string.
    """
    iterate_module = importlib.import_module("chemsmart.cli.iterate.iterate")

    def _patch(run_return="out.xyz"):
        mock_job = MagicMock()
        mock_job.run.return_value = run_return
        mocker.patch.object(
            iterate_module, "IterateJob", return_value=mock_job
        )
        mocker.patch.object(iterate_module, "IterateJobRunner")
        mocker.patch.object(iterate_module, "IterateJobSettings")
        mocker.patch.object(
            iterate_module,
            "validate_config",
            return_value={"skeletons": [1], "substituents": [1]},
        )
        return mock_job

    return _patch
