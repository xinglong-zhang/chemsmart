"""Focused tests for the PySCF Click surface and option propagation."""

from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from chemsmart.cli.pyscf import pyscf
from chemsmart.cli.pyscf.common import resolve_engine, settings_for_molecule
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from tests.helpers import run_pyscf_and_capture_settings


@pytest.mark.parametrize("jobtype", ["sp", "opt", "hess"])
def test_group_and_all_leaves_expose_help(single_molecule_xyz_file, jobtype):
    runner = CliRunner()

    group = runner.invoke(pyscf, ["--help"])
    leaf = runner.invoke(
        pyscf,
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            jobtype,
            "--help",
        ],
    )

    assert group.exit_code == 0, group.output
    assert leaf.exit_code == 0, leaf.output
    assert jobtype in group.output


@pytest.mark.parametrize(
    ("option", "field", "expected"),
    [
        (["-t", "receipt"], "title", "receipt"),
        (["-c", "1"], "charge", 1),
        (["-m", "2"], "multiplicity", 2),
        (["-x", "pbe0"], "functional", "pbe0"),
        (["-A", "hf"], "ab_initio", "hf"),
        (["-b", "def2-tzvp"], "basis", "def2-tzvp"),
        (
            ["--aux-basis", "def2-universal-jkfit"],
            "aux_basis",
            "def2-universal-jkfit",
        ),
        (["-d", "d3bj"], "dispersion", "d3bj"),
        (["--defgrid", "defgrid3"], "defgrid", "defgrid3"),
        (["--scf-tol", "1e-10"], "scf_tol", 1e-10),
        (["--scf-maxiter", "75"], "scf_maxiter", 75),
        (["--no-density-fit"], "density_fit", False),
        (["--opt-solver", "geometric"], "opt_solver", "geometric"),
        (["--opt-maxsteps", "40"], "opt_maxsteps", 40),
        (["-sm", "cpcm"], "solvent_model", "cpcm"),
        (["-si", "toluene"], "solvent_id", "toluene"),
    ],
)
def test_every_group_override_survives_the_merge(
    single_molecule_xyz_file, option, field, expected
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            *option,
            "sp",
        ]
    )

    assert result.exit_code == 0, result.output
    assert settings is not None
    assert getattr(settings, field) == expected

    if field == "functional":
        assert settings.ab_initio is None
    if field == "ab_initio":
        assert settings.functional is None


@pytest.mark.parametrize(
    ("option", "expected"), [(["--gpu"], "gpu"), (["--no-gpu"], "cpu")]
)
def test_explicit_engine_override_wins(
    single_molecule_xyz_file, option, expected
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            *option,
            "sp",
        ],
        num_gpus=1 if expected == "cpu" else 0,
    )

    assert result.exit_code == 0, result.output
    assert settings.engine == expected


@pytest.mark.parametrize(("num_gpus", "expected"), [(0, "cpu"), (1, "gpu")])
def test_server_gpu_count_resolves_default_engine(
    single_molecule_xyz_file, num_gpus, expected
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "sp",
        ],
        num_gpus=num_gpus,
    )

    assert result.exit_code == 0, result.output
    assert settings.engine == expected


def test_remove_solvent_clears_project_single_point_solvent(
    single_molecule_xyz_file,
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "--remove-solvent",
            "sp",
        ]
    )

    assert result.exit_code == 0, result.output
    assert (settings.solvent_model, settings.solvent_id) == (None, None)


@pytest.mark.parametrize("jobtype", ["sp", "opt", "hess"])
def test_leaf_selects_matching_project_settings(
    single_molecule_xyz_file, jobtype
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            jobtype,
        ],
        jobtype=jobtype,
    )

    assert result.exit_code == 0, result.output
    assert settings.jobtype == jobtype
    assert settings.project_yaml_digest is not None


@pytest.mark.parametrize(
    ("project_engine", "num_gpus", "expected"),
    [
        ("gpu", None, "gpu"),
        ("cpu", 2, "cpu"),
        (None, 2, "gpu"),
        (None, 0, "cpu"),
    ],
)
def test_project_engine_is_respected_before_server_gpus(
    project_engine, num_gpus, expected
):
    settings = PySCFJobSettings(
        functional="b3lyp", basis="def2-svp", engine=project_engine
    )
    ctx = SimpleNamespace(
        obj={"jobrunner": SimpleNamespace(num_gpus=num_gpus)}
    )
    assert resolve_engine(ctx, settings) == expected


def test_charge_and_multiplicity_fall_back_to_the_molecule():
    settings = PySCFJobSettings(functional="b3lyp", basis="def2-svp")
    molecule = SimpleNamespace(charge=1, multiplicity=2)
    resolved = settings_for_molecule(settings, molecule)
    assert (resolved.charge, resolved.multiplicity) == (1, 2)
    assert settings.charge is None


def test_explicit_charge_and_multiplicity_win_over_the_molecule():
    settings = PySCFJobSettings(
        functional="b3lyp", basis="def2-svp", charge=0, multiplicity=1
    )
    molecule = SimpleNamespace(charge=1, multiplicity=2)
    resolved = settings_for_molecule(settings, molecule)
    assert (resolved.charge, resolved.multiplicity) == (0, 1)


def test_missing_charge_and_multiplicity_are_reported():
    settings = PySCFJobSettings(functional="b3lyp", basis="def2-svp")
    molecule = SimpleNamespace(charge=None, multiplicity=None)
    with pytest.raises(ValueError, match="Charge and multiplicity"):
        settings_for_molecule(settings, molecule)
