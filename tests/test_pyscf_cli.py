"""Focused tests for the PySCF Click surface and option propagation."""

import glob
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.pyscf import pyscf
from chemsmart.cli.pyscf.common import resolve_engine, settings_for_molecule
from chemsmart.cli.run import run
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import read_pyscf_h5
from chemsmart.jobs.pyscf.runner import PySCFJobRunner
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from tests.helpers import run_pyscf_and_capture_settings


@pytest.mark.parametrize("jobtype", ["sp", "opt", "hess"])
def test_group_and_all_leaves_expose_help(pyscf_water_xyz, jobtype):
    runner = CliRunner()

    group = runner.invoke(pyscf, ["--help"])
    leaf = runner.invoke(
        pyscf,
        [
            "-p",
            "test",
            "-f",
            pyscf_water_xyz,
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
        (["-c", "2"], "charge", 2),
        (["-m", "3"], "multiplicity", 3),
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
    pyscf_water_xyz, option, field, expected
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            pyscf_water_xyz,
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
def test_explicit_engine_override_wins(pyscf_water_xyz, option, expected):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            pyscf_water_xyz,
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
    pyscf_water_xyz, num_gpus, expected
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            pyscf_water_xyz,
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
    pyscf_water_xyz,
):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            pyscf_water_xyz,
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
def test_leaf_selects_matching_project_settings(pyscf_water_xyz, jobtype):
    result, settings = run_pyscf_and_capture_settings(
        [
            "-p",
            "test",
            "-f",
            pyscf_water_xyz,
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
    molecule = Molecule(
        symbols=["He"], positions=np.zeros((1, 3)), charge=1, multiplicity=2
    )
    resolved = settings_for_molecule(settings, molecule)
    assert (resolved.charge, resolved.multiplicity) == (1, 2)
    assert settings.charge is None


def test_explicit_charge_and_multiplicity_win_over_the_molecule():
    settings = PySCFJobSettings(
        functional="b3lyp", basis="def2-svp", charge=0, multiplicity=1
    )
    molecule = Molecule(
        symbols=["He"], positions=np.zeros((1, 3)), charge=1, multiplicity=2
    )
    resolved = settings_for_molecule(settings, molecule)
    assert (resolved.charge, resolved.multiplicity) == (0, 1)


def test_missing_charge_and_multiplicity_are_reported():
    settings = PySCFJobSettings(functional="b3lyp", basis="def2-svp")
    molecule = Molecule(symbols=["He"], positions=np.zeros((1, 3)))
    with pytest.raises(ValueError, match="Charge and multiplicity"):
        settings_for_molecule(settings, molecule)


class TestPySCFCliRequests:
    def _invoke(self, args):
        return CliRunner().invoke(
            pyscf,
            args,
            obj={"jobrunner": MagicMock(num_gpus=0)},
            catch_exceptions=False,
        )

    @pytest.mark.parametrize(
        ("options", "message"),
        [
            (["-x", "pbe", "-A", "hf"], "either --functional or --ab-initio"),
            (["-A", "hf", "--defgrid", "defgrid2"], "DFT-only"),
            (["--no-density-fit", "-ab", "x"], "--aux-basis cannot"),
        ],
    )
    def test_contradictory_options_are_refused(
        self, pyscf_water_xyz, options, message
    ):
        with pytest.raises(ValueError, match=message):
            self._invoke(["-p", "test", "-f", pyscf_water_xyz, *options, "sp"])

    def test_a_structure_source_is_required(self):
        with pytest.raises(ValueError, match="has not been specified"):
            self._invoke(["-p", "test", "sp"])

    def test_only_one_structure_source_is_allowed(self, pyscf_water_xyz):
        with pytest.raises(ValueError, match="Both"):
            self._invoke(
                ["-p", "test", "-f", pyscf_water_xyz, "-P", "962", "sp"]
            )

    def test_label_and_append_label_conflict(self, pyscf_water_xyz):
        with pytest.raises(ValueError, match="but not both"):
            self._invoke(
                [
                    "-p",
                    "test",
                    "-f",
                    pyscf_water_xyz,
                    "-l",
                    "a",
                    "-a",
                    "b",
                    "sp",
                ]
            )


class TestPySCFCliRealRuns:
    """Preview runs through ``chemsmart run`` from real input files."""

    def _preview(self, *args):
        result = CliRunner().invoke(
            run, ["--fake", "pyscf", "-p", "test", *args], obj={}
        )
        assert result.exit_code == 0, (result.output, result.exception)

    def test_append_label_extends_the_file_name(
        self, pyscf_water_xyz, temporary_working_dir
    ):
        self._preview(
            "-f", pyscf_water_xyz, "-c", "0", "-m", "1", "-a", "extra", "sp"
        )
        assert glob.glob("water_extra*.h5")

    def test_forces_and_tolerance_reach_the_results(
        self, pyscf_water_xyz, temporary_working_dir
    ):
        self._preview(
            "-f",
            pyscf_water_xyz,
            "-c",
            "0",
            "-m",
            "1",
            "--forces",
            "--point-group-tol",
            "0.01",
            "-l",
            "requested",
            "sp",
        )
        spec = read_pyscf_h5(glob.glob("requested*.h5")[0])[0]
        assert spec["forces"] is True
        assert spec["point_group_tol"] == 0.01

    def test_no_density_fit_clears_an_inherited_auxiliary_basis(
        self, pyscf_water_xyz, temporary_working_dir
    ):
        self._preview(
            "-f",
            pyscf_water_xyz,
            "-c",
            "0",
            "-m",
            "1",
            "--no-density-fit",
            "-l",
            "nodf",
            "sp",
        )
        spec = read_pyscf_h5(glob.glob("nodf*.h5")[0])[0]
        assert spec["density_fit"] is False
        assert spec["aux_basis"] is None

    def test_each_selected_structure_gets_its_own_job(
        self, multiple_molecules_xyz_file, temporary_working_dir
    ):
        self._preview(
            "-f",
            multiple_molecules_xyz_file,
            "-c",
            "0",
            "-m",
            "2",
            "-i",
            "1,2",
            "sp",
        )
        assert len(glob.glob("*_idx1_*.h5")) == 1
        assert len(glob.glob("*_idx2_*.h5")) == 1

    def test_single_selected_index_runs_one_job(
        self, multiple_molecules_xyz_file, temporary_working_dir
    ):
        self._preview(
            "-f",
            multiple_molecules_xyz_file,
            "-c",
            "0",
            "-m",
            "2",
            "-i",
            "2",
            "sp",
        )
        assert len(glob.glob("*.h5")) == 1

    def test_append_label_works_for_a_pubchem_structure(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        with patch.object(
            Molecule, "from_pubchem", return_value=[pyscf_water_molecule]
        ):
            self._preview(
                "-P", "962", "-c", "0", "-m", "1", "-a", "extra", "sp"
            )
        assert glob.glob("output_extra*.h5")

    def test_pubchem_structure_is_used_with_a_default_label(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        with patch.object(
            Molecule, "from_pubchem", return_value=[pyscf_water_molecule]
        ):
            self._preview("-P", "962", "-c", "0", "-m", "1", "sp")
        assert glob.glob("output_sp*.h5")


class TestImpossibleElectronStates:
    """Electron counts PySCF would refuse, stopped before any job starts."""

    def test_odd_electron_count_cannot_be_a_singlet(
        self, pyscf_oh_molecule, temporary_working_dir
    ):
        settings = PySCFJobSettings(
            functional="b3lyp", basis="def2-svp", charge=0, multiplicity=1
        )
        with pytest.raises(ValueError, match="needs an even multiplicity"):
            settings_for_molecule(settings, pyscf_oh_molecule)

    def test_the_message_says_what_to_check(self, pyscf_oh_molecule):
        settings = PySCFJobSettings(
            functional="b3lyp", basis="def2-svp", charge=0, multiplicity=1
        )
        with pytest.raises(ValueError, match="9 electrons.*Check -c"):
            settings_for_molecule(settings, pyscf_oh_molecule)

    @pytest.mark.parametrize(
        ("symbols", "charge", "multiplicity", "message"),
        [
            (["H", "H"], 0, 2, "needs an odd multiplicity"),
            (["H", "H"], 0, 5, "cannot have 4 unpaired"),
            (["He"], 3, 2, "negative number of electrons"),
        ],
    )
    def test_impossible_states(self, symbols, charge, multiplicity, message):
        molecule = Molecule(
            symbols=symbols,
            positions=np.arange(len(symbols) * 3.0).reshape(-1, 3),
        )
        settings = PySCFJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            charge=charge,
            multiplicity=multiplicity,
        )
        with pytest.raises(ValueError, match=message):
            settings.validate_electrons(molecule)

    @pytest.mark.parametrize(
        ("symbols", "charge", "multiplicity"),
        [
            (["O", "H"], 0, 2),
            (["O", "H"], -1, 1),
            (["H", "H"], 0, 3),
            (["H"], 1, 1),
        ],
    )
    def test_possible_states_are_accepted(self, symbols, charge, multiplicity):
        molecule = Molecule(
            symbols=symbols,
            positions=np.arange(len(symbols) * 3.0).reshape(-1, 3),
        )
        settings = PySCFJobSettings(
            functional="b3lyp",
            basis="sto-3g",
            charge=charge,
            multiplicity=multiplicity,
        )
        settings.validate_electrons(molecule)

    def test_fake_run_stops_before_writing_anything(
        self, pyscf_oh_molecule, temporary_working_dir
    ):
        with patch.object(
            Molecule, "from_pubchem", return_value=[pyscf_oh_molecule]
        ):
            result = CliRunner().invoke(
                run,
                ["--fake", "pyscf", "-p", "test", "-P", "962"]
                + ["-c", "0", "-m", "1", "sp"],
                obj={},
            )
        assert result.exit_code != 0
        assert "needs an even multiplicity" in str(result.exception)
        assert glob.glob("*.py") == []


class TestFailedJobsThroughTheCli:
    def _run(self, *args):
        return CliRunner().invoke(run, ["pyscf", "-p", "test", *args], obj={})

    def test_real_failure_exits_non_zero_and_gives_the_reason(
        self, pyscf_water_xyz, temporary_working_dir
    ):
        pytest.importorskip("pyscf")
        result = self._run(
            "-f",
            pyscf_water_xyz,
            "-c",
            "0",
            "-m",
            "1",
            "-si",
            "not_a_solvent",
            "-l",
            "bad",
            "sp",
        )
        assert result.exit_code == 1
        assert "not_a_solvent" in result.output
        assert "failed" in result.output

    def test_batch_runs_every_structure_before_reporting(
        self, multiple_molecules_xyz_file, temporary_working_dir
    ):
        dies = f"{sys.executable} -c 'raise SystemExit(1)'"
        with patch.object(PySCFJobRunner, "_get_command", return_value=dies):
            result = self._run(
                "-f",
                multiple_molecules_xyz_file,
                "-c",
                "0",
                "-m",
                "2",
                "-i",
                "1,2",
                "sp",
            )
        assert result.exit_code == 1
        assert "2 of 2 jobs failed" in result.output
        assert len(glob.glob("*_idx1_*.py")) == 1
        assert len(glob.glob("*_idx2_*.py")) == 1

    def test_passing_preview_still_exits_zero(
        self, pyscf_water_xyz, temporary_working_dir
    ):
        result = CliRunner().invoke(
            run,
            ["--fake", "pyscf", "-p", "test", "-f", pyscf_water_xyz]
            + ["-c", "0", "-m", "1", "sp"],
            obj={},
        )
        assert result.exit_code == 0, (result.output, result.exception)
