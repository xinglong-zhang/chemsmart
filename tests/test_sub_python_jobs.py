"""``chemsmart sub`` for jobs that run in Python, with no external program.

These commands calculate as they parse under ``chemsmart run``. Under
``chemsmart sub`` they must only write scripts, and those scripts must repeat
the same calculation on the compute node.
"""

import os
import shutil

import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from tests.helpers import (
    assert_uses_submitting_python,
    read_grouper_structures,
    run_generated_script_command,
    submit_for_test,
)


class TestThermochemistrySubmission:
    def test_writes_scripts_without_calculating(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        gaussian_co2_opt_outfile,
    ):
        monkeypatch.chdir(tmp_path)
        result = submit_for_test(
            server_yaml_file,
            "thermochemistry",
            "-T",
            298.15,
            "-f",
            gaussian_co2_opt_outfile,
        )
        assert result.exit_code == 0, result.output

        label = os.path.splitext(os.path.basename(gaussian_co2_opt_outfile))[0]
        assert (
            tmp_path / f"chemsmart_run_thermochemistry_{label}.py"
        ).exists()
        assert_uses_submitting_python(
            tmp_path / f"chemsmart_sub_thermochemistry_{label}.sh"
        )
        assert not list(tmp_path.glob("*.dat"))

    def test_one_script_covers_all_files_and_matches_a_local_run(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        gaussian_co2_opt_outfile,
        gaussian_conformer1_outfile,
    ):
        """Several files are one job, with options kept on the node."""
        args = [
            "thermochemistry",
            "-T",
            "298.15",
            "-csg",
            "100",
            "-f",
            gaussian_co2_opt_outfile,
            "-f",
            gaussian_conformer1_outfile,
            "-o",
            "thermo.dat",
        ]
        local = tmp_path / "local"
        local.mkdir()
        monkeypatch.chdir(local)
        assert CliRunner().invoke(run, args).exit_code == 0

        node = tmp_path / "node"
        node.mkdir()
        monkeypatch.chdir(node)
        result = submit_for_test(server_yaml_file, *args)
        assert result.exit_code == 0, result.output
        run_scripts = list(node.glob("chemsmart_run_*.py"))
        assert len(run_scripts) == 1

        node_result = run_generated_script_command(run_scripts[0])
        assert node_result.exit_code == 0, node_result.output
        assert (node / "thermo.dat").read_text() == (
            local / "thermo.dat"
        ).read_text()

    def test_no_files_is_an_error(
        self, tmp_path, monkeypatch, server_yaml_file
    ):
        monkeypatch.chdir(tmp_path)
        result = submit_for_test(
            server_yaml_file, "thermochemistry", "-T", 298.15
        )
        assert result.exit_code == 1
        assert "No thermochemistry output files to submit" in result.output
        assert not list(tmp_path.glob("chemsmart_*"))


class TestBoltzmannSubmission:
    def test_submits_and_then_gives_the_local_result(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        gaussian_conformer1_outfile,
        gaussian_conformer2_outfile,
    ):
        args = [
            "thermochemistry",
            "-T",
            "298.15",
            "-f",
            gaussian_conformer1_outfile,
            "-f",
            gaussian_conformer2_outfile,
            "boltzmann",
            "-w",
            "electronic",
        ]
        local = tmp_path / "local"
        local.mkdir()
        monkeypatch.chdir(local)
        assert CliRunner().invoke(run, args).exit_code == 0
        local_reports = list(local.glob("*.dat"))
        assert len(local_reports) == 1

        node = tmp_path / "node"
        node.mkdir()
        monkeypatch.chdir(node)
        result = submit_for_test(server_yaml_file, *args)
        assert result.exit_code == 0, result.output
        assert not list(node.glob("*.dat"))
        run_scripts = list(node.glob("chemsmart_run_*.py"))
        assert len(run_scripts) == 1
        assert "'boltzmann'" in run_scripts[0].read_text()

        node_result = run_generated_script_command(run_scripts[0])
        assert node_result.exit_code == 0, node_result.output
        node_report = node / local_reports[0].name
        assert node_report.read_text() == local_reports[0].read_text()


class TestIterateSubmission:
    # Few samples keep the run short; these flags must also survive rebuilding
    SAMPLING = ["-s", "16", "-a", "2"]

    def test_writes_scripts_without_generating_structures(
        self, iterate_submission_folder, server_yaml_file
    ):
        folder = iterate_submission_folder
        result = submit_for_test(
            server_yaml_file, "iterate", "-f", "config.toml", *self.SAMPLING
        )
        assert result.exit_code == 0, result.output
        assert (folder / "chemsmart_run_iterate_config.py").exists()
        assert_uses_submitting_python(
            folder / "chemsmart_sub_iterate_config.sh"
        )
        assert not (folder / "iterate_out.xyz").exists()

    def test_generated_command_matches_a_local_run(
        self, iterate_submission_folder, server_yaml_file
    ):
        folder = iterate_submission_folder
        args = ["iterate", "-f", "config.toml", *self.SAMPLING]
        assert (
            CliRunner().invoke(run, [*args, "-o", "local"]).exit_code == 0
        ), "local run failed"

        result = submit_for_test(server_yaml_file, *args, "-o", "node")
        assert result.exit_code == 0, result.output
        node_result = run_generated_script_command(
            folder / "chemsmart_run_iterate_config.py"
        )
        assert node_result.exit_code == 0, node_result.output
        assert (folder / "node.xyz").read_text() == (
            folder / "local.xyz"
        ).read_text()

    def test_separate_outputs_survive_the_generated_command(
        self, iterate_submission_folder, server_yaml_file
    ):
        """The default output name must not reach the node as an option."""
        folder = iterate_submission_folder
        (folder / "separate").mkdir()
        result = submit_for_test(
            server_yaml_file,
            "iterate",
            "-f",
            "config.toml",
            "--separate-outputs",
            "-d",
            "separate",
            *self.SAMPLING,
        )
        assert result.exit_code == 0, result.output

        node_result = run_generated_script_command(
            folder / "chemsmart_run_iterate_config.py"
        )
        assert node_result.exit_code == 0, node_result.output
        assert list((folder / "separate").glob("*.xyz"))


class TestGrouperSubmission:
    def test_writes_scripts_without_grouping(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        multiple_molecules_xyz_file,
    ):
        monkeypatch.chdir(tmp_path)
        result = submit_for_test(
            server_yaml_file,
            "grouper",
            "-f",
            multiple_molecules_xyz_file,
            "rmsd",
        )
        assert result.exit_code == 0, result.output
        assert (
            tmp_path / "chemsmart_run_grouper_crest_conformers_rmsd.py"
        ).exists()
        assert_uses_submitting_python(
            tmp_path / "chemsmart_sub_grouper_crest_conformers_rmsd.sh"
        )
        assert not list(tmp_path.glob("*_group_result"))

    def test_generated_command_matches_a_local_run(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        multiple_molecules_xyz_file,
    ):
        """Every non-default option changes the result or its file names."""
        args = [
            "grouper",
            "-f",
            multiple_molecules_xyz_file,
            "-ih",
            "-np",
            "2",
            "-N",
            "3",
            "-m",
            "csv",
            "-r",
            "center",
            "-a",
            "tagged",
            "rmsd",
        ]
        local = tmp_path / "local"
        local.mkdir()
        monkeypatch.chdir(local)
        assert CliRunner().invoke(run, args).exit_code == 0
        local_structures = read_grouper_structures(local)
        assert len(local_structures) == 3

        node = tmp_path / "node"
        node.mkdir()
        monkeypatch.chdir(node)
        result = submit_for_test(server_yaml_file, *args)
        assert result.exit_code == 0, result.output
        run_scripts = list(node.glob("chemsmart_run_*.py"))
        assert len(run_scripts) == 1

        node_result = run_generated_script_command(run_scripts[0])
        assert node_result.exit_code == 0, node_result.output
        assert read_grouper_structures(node) == local_structures
        assert list(node.rglob("*.csv"))

    def test_directory_mode_keeps_thermochemistry_options(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        ts_conformers_log_directory,
    ):
        """The energies written depend on the options reaching the node."""
        args = [
            "grouper",
            "-d",
            ".",
            "-p",
            "gaussian",
            "-E",
            "qhG",
            "-temp",
            "350",
            "-c",
            "0.5",
            "--no-weighted",
            "--no-check-imaginary-frequencies",
            "rmsd",
        ]
        local = tmp_path / "local"
        shutil.copytree(ts_conformers_log_directory, local)
        monkeypatch.chdir(local)
        assert CliRunner().invoke(run, args).exit_code == 0
        local_structures = read_grouper_structures(local)
        assert local_structures
        assert all("qhG(Hartree)" in s for s in local_structures.values())

        node = tmp_path / "node"
        shutil.copytree(ts_conformers_log_directory, node)
        monkeypatch.chdir(node)
        result = submit_for_test(server_yaml_file, *args)
        assert result.exit_code == 0, result.output
        assert not read_grouper_structures(node)
        run_scripts = list(node.glob("chemsmart_run_*.py"))
        assert len(run_scripts) == 1

        node_result = run_generated_script_command(run_scripts[0])
        assert node_result.exit_code == 0, node_result.output
        assert read_grouper_structures(node) == local_structures

    @pytest.mark.parametrize(
        "strategy_args, expected",
        [
            (["irmsd", "--inversion", "off"], "'--inversion', 'off'"),
            (
                ["tanimoto", "-ft", "morgan"],
                "'--fingerprint-type', 'morgan'",
            ),
            (["tfd", "--use-weights"], "'--use-weights'"),
            (["pymolrmsd"], "'pymolrmsd'"),
        ],
    )
    def test_strategy_options_reach_the_node(
        self,
        strategy_args,
        expected,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        multiple_molecules_xyz_file,
    ):
        """Checks the command only, as some strategies need optional tools."""
        monkeypatch.chdir(tmp_path)
        result = submit_for_test(
            server_yaml_file,
            "grouper",
            "-f",
            multiple_molecules_xyz_file,
            *strategy_args,
        )
        assert result.exit_code == 0, result.output
        run_scripts = list(tmp_path.glob("chemsmart_run_*.py"))
        assert len(run_scripts) == 1
        assert expected in run_scripts[0].read_text()
