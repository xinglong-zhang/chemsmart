"""Tests for ``chemsmart run sterimol`` and ``chemsmart sub sterimol``.

Reference values come from morfeus-ml 0.8.0 (Sterimol with RDKit radii, 36000
rotation vectors, L_value_uncorrected, and the same excluded atoms). L and B5
are analytical. B1 is sampled more coarsely here, so it has a wider tolerance.
"""

import ast
import json
import os
import re
import sys

import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.sterimol.runner import SterimolJobRunner
from tests.helpers import (
    assert_sterimol_matches,
    invoke_sterimol,
    run_sterimol_to_json,
)

MORFEUS_REFERENCE = {
    (2, 4): (5.75827, 1.76204, 3.85632),
    (4, 2): (4.43546, 2.11960, 3.38237),
    (4, 5): (3.93579, 1.88007, 3.47812),
    (2, 3): (3.23416, 1.70809, 2.24541),
}


class TestSterimolValues:
    @pytest.mark.parametrize("axis, reference", MORFEUS_REFERENCE.items())
    def test_auto_fragment_matches_morfeus(
        self, tmp_path, sterimol_halo_chain_xyz_file, axis, reference
    ):
        report = run_sterimol_to_json(
            tmp_path,
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            *axis,
            "--auto-fragment",
        )
        assert_sterimol_matches(report, reference)

    def test_explicit_atoms_match_auto_fragment(
        self, tmp_path, sterimol_halo_chain_xyz_file
    ):
        """Ranges are 1-based, and order and spacing do not matter."""
        report = run_sterimol_to_json(
            tmp_path,
            "--file",
            sterimol_halo_chain_xyz_file,
            "--axis",
            2,
            4,
            "--atoms",
            "15-17, 4-7",
        )
        assert report["selection"] == "explicit"
        assert report["substituent_atoms"] == [4, 5, 6, 7, 15, 16, 17]
        assert_sterimol_matches(report, MORFEUS_REFERENCE[(2, 4)])

    def test_explicit_atoms_work_for_a_ring_bond(
        self, tmp_path, sterimol_toluene_xyz_file
    ):
        """--atoms is the way to measure what auto-fragment refuses."""
        report = run_sterimol_to_json(
            tmp_path,
            "-f",
            sterimol_toluene_xyz_file,
            "--axis",
            2,
            3,
            "--atoms",
            "3-5,12,13",
        )
        assert report["substituent_atoms"] == [3, 4, 5, 12, 13]
        assert report["L"] > 0


class TestSterimolIndex:
    def test_default_is_the_last_structure(
        self, tmp_path, sterimol_trajectory_xyz_file
    ):
        args = [
            "-f",
            sterimol_trajectory_xyz_file,
            "--axis",
            2,
            4,
            "--auto-fragment",
        ]
        first = run_sterimol_to_json(tmp_path, *args, "-i", 1)
        second = run_sterimol_to_json(tmp_path, *args, "-i", 2)
        default = run_sterimol_to_json(tmp_path, *args)

        assert_sterimol_matches(first, MORFEUS_REFERENCE[(2, 4)])
        assert default["B5"] == pytest.approx(second["B5"])
        assert default["B5"] != pytest.approx(first["B5"], abs=1e-3)

    def test_rejects_several_structures(self, sterimol_trajectory_xyz_file):
        result = invoke_sterimol(
            "-f",
            sterimol_trajectory_xyz_file,
            "--axis",
            2,
            4,
            "--auto-fragment",
            "-i",
            ":",
        )
        assert result.exit_code == 1
        assert "Select a single structure" in result.output

    def test_rejects_missing_structure(self, sterimol_trajectory_xyz_file):
        result = invoke_sterimol(
            "-f",
            sterimol_trajectory_xyz_file,
            "--axis",
            2,
            4,
            "--auto-fragment",
            "-i",
            5,
        )
        assert result.exit_code == 1
        assert "--index 5 does not select a structure" in result.output


class TestSterimolReport:
    def test_json_records_how_the_values_were_obtained(
        self, tmp_path, sterimol_halo_chain_xyz_file
    ):
        report = run_sterimol_to_json(
            tmp_path,
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            4,
            5,
            "--auto-fragment",
        )
        assert report["file"] == sterimol_halo_chain_xyz_file
        assert report["index"] == "-1"
        assert report["axis"] == {"atoms": [4, 5], "symbols": ["C", "C"]}
        assert report["selection"] == "auto-fragment"
        assert report["substituent_atoms"] == [5, 6, 7, 17]
        assert report["units"] == "Angstrom"
        assert "RDKit" in report["radius_source"]
        assert "+0.40" in report["length_convention"]

    def test_terminal_output(
        self, tmp_path, monkeypatch, sterimol_halo_chain_xyz_file
    ):
        monkeypatch.chdir(tmp_path)
        result = invoke_sterimol(
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            2,
            4,
            "--auto-fragment",
        )
        assert result.exit_code == 0, result.output
        assert "Axis:         C2 -> C4" in result.output
        assert "Substituent:  4-7,15-17 (auto-fragment)" in result.output
        assert "L  = 5.758 Å" in result.output
        assert "B1 = 1.763 Å" in result.output
        assert "B5 = 3.856 Å" in result.output


class TestSterimolErrors:
    @pytest.mark.parametrize(
        "selection",
        [[], ["--auto-fragment", "--atoms", "4"]],
        ids=["neither", "both"],
    )
    def test_needs_exactly_one_selection(
        self, sterimol_halo_chain_xyz_file, selection
    ):
        result = invoke_sterimol(
            "-f", sterimol_halo_chain_xyz_file, "--axis", 2, 4, *selection
        )
        assert result.exit_code == 2
        assert "exactly one of --auto-fragment or --atoms" in result.output

    @pytest.mark.parametrize(
        "axis, message",
        [
            ((2, 3), "is in a ring"),
            ((1, 4), "are not bonded"),
        ],
    )
    def test_auto_fragment_points_to_atoms_option(
        self, sterimol_toluene_xyz_file, axis, message
    ):
        result = invoke_sterimol(
            "-f",
            sterimol_toluene_xyz_file,
            "--axis",
            *axis,
            "--auto-fragment",
        )
        assert result.exit_code == 1
        assert message in result.output
        assert "Select the atoms with --atoms" in result.output

    @pytest.mark.parametrize(
        "axis, atoms, message",
        [
            ((2, 4), "2,4", "cannot contain the first axis atom"),
            ((2, 4), "5,6", "must contain the second axis atom"),
            ((2, 4), "4,99", r"\[99\] are outside"),
            ((2, 99), "4", r"\[99\] are outside"),
            ((2, 4), "abc", "abc"),
        ],
    )
    def test_rejects_bad_explicit_atoms(
        self, sterimol_halo_chain_xyz_file, axis, atoms, message
    ):
        result = invoke_sterimol(
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            *axis,
            "--atoms",
            atoms,
        )
        assert result.exit_code == 1
        assert re.search(message, result.output)

    def test_rejects_zero_based_axis(self, sterimol_halo_chain_xyz_file):
        """Index 0 is refused before any file is read."""
        result = invoke_sterimol(
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            0,
            1,
            "--auto-fragment",
        )
        assert result.exit_code == 2

    def test_auto_fragment_out_of_range_has_no_atoms_hint(
        self, sterimol_halo_chain_xyz_file
    ):
        """--atoms would not help with a wrong axis, so it is not suggested."""
        result = invoke_sterimol(
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            2,
            99,
            "--auto-fragment",
        )
        assert result.exit_code == 1
        assert "outside the valid 1-based range" in result.output
        assert "--atoms" not in result.output

    def test_rejects_an_empty_file(self, tmp_path):
        empty_file = tmp_path / "empty.xyz"
        empty_file.touch()
        result = invoke_sterimol(
            "-f", empty_file, "--axis", 1, 2, "--auto-fragment"
        )
        assert result.exit_code == 1
        assert "No structure could be read" in result.output


class TestSterimolRunner:
    def test_needs_no_executable(self, pbs_server):
        """It must not look for a program section in server.yaml."""
        assert SterimolJobRunner(server=pbs_server).executable is None


class TestSterimolOutputFile:
    def test_default_name_uses_file_and_axis(
        self, tmp_path, monkeypatch, sterimol_halo_chain_xyz_file
    ):
        monkeypatch.chdir(tmp_path)
        args = ["-f", sterimol_halo_chain_xyz_file, "--auto-fragment"]
        assert (
            invoke_sterimol(*args, "--axis", 2, 4).exit_code == 0
            and invoke_sterimol(*args, "--axis", 4, 2).exit_code == 0
        )

        report = json.loads((tmp_path / "chain_2_4_sterimol.json").read_text())
        assert_sterimol_matches(report, MORFEUS_REFERENCE[(2, 4)])
        report = json.loads((tmp_path / "chain_4_2_sterimol.json").read_text())
        assert_sterimol_matches(report, MORFEUS_REFERENCE[(4, 2)])

    def test_default_name_uses_selected_structure(
        self, tmp_path, monkeypatch, sterimol_trajectory_xyz_file
    ):
        monkeypatch.chdir(tmp_path)
        result = invoke_sterimol(
            "-f",
            sterimol_trajectory_xyz_file,
            "--axis",
            2,
            4,
            "--auto-fragment",
            "-i",
            1,
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "traj_2_4_i1_sterimol.json").exists()

    def test_label_names_the_report(
        self, tmp_path, monkeypatch, sterimol_halo_chain_xyz_file
    ):
        monkeypatch.chdir(tmp_path)
        result = invoke_sterimol(
            "-f",
            sterimol_halo_chain_xyz_file,
            "--axis",
            2,
            4,
            "--auto-fragment",
            "-l",
            "isobutyl",
        )
        assert result.exit_code == 0, result.output
        assert (tmp_path / "isobutyl_sterimol.json").exists()

    def test_rerun_replaces_an_existing_report(
        self, tmp_path, sterimol_halo_chain_xyz_file
    ):
        """The same axis with other atoms must not reuse the old values."""
        args = ["-f", sterimol_halo_chain_xyz_file, "--axis", 2, 4, "--atoms"]
        first = run_sterimol_to_json(tmp_path, *args, "4-7,15-17")
        second = run_sterimol_to_json(tmp_path, *args, "4,15,16")

        assert first["substituent_atoms"] == [4, 5, 6, 7, 15, 16, 17]
        assert second["substituent_atoms"] == [4, 15, 16]
        assert second["L"] != pytest.approx(first["L"])


class TestSterimolSubmission:
    """``sub`` writes scripts that repeat the command on a compute node."""

    def test_writes_scripts_without_a_program_section(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        sterimol_halo_chain_xyz_file,
    ):
        """The node finds the submitting Python, not a configured program."""
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "sterimol",
                "-f",
                sterimol_halo_chain_xyz_file,
                "--axis",
                "2",
                "4",
                "--auto-fragment",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output

        run_script = tmp_path / "chemsmart_run_sterimol_chain_2_4.py"
        submit_script = tmp_path / "chemsmart_sub_sterimol_chain_2_4.sh"
        assert "'--test'" not in run_script.read_text()
        submit_text = submit_script.read_text()
        assert f"export PATH={os.path.dirname(sys.executable)}:" in submit_text
        assert "#PBS -N sterimol_chain_2_4" in submit_text
        assert "conda activate" not in submit_text
        # submitting calculates nothing
        assert not (tmp_path / "chain_2_4_sterimol.json").exists()

    @pytest.mark.parametrize(
        "selection",
        [["--auto-fragment"], ["--atoms", "4-7,15-17"]],
        ids=["auto-fragment", "explicit"],
    )
    def test_generated_command_gives_the_local_values(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        sterimol_halo_chain_xyz_file,
        selection,
    ):
        """Every option must survive the rebuilt command line."""
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "sterimol",
                "-f",
                sterimol_halo_chain_xyz_file,
                "--axis",
                "2",
                "4",
                "-i",
                "1",
                *selection,
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output

        run_script = tmp_path / "chemsmart_run_sterimol_chain_2_4_i1.py"
        match = re.search(r"run\((\[.*\])\)", run_script.read_text())
        node_result = CliRunner().invoke(run, ast.literal_eval(match.group(1)))
        assert node_result.exit_code == 0, node_result.output

        report = json.loads(
            (tmp_path / "chain_2_4_i1_sterimol.json").read_text()
        )
        assert report["index"] == "1"
        assert report["substituent_atoms"] == [4, 5, 6, 7, 15, 16, 17]
        assert_sterimol_matches(report, MORFEUS_REFERENCE[(2, 4)])

    def test_bad_selection_fails_before_submission(
        self,
        tmp_path,
        monkeypatch,
        server_yaml_file,
        sterimol_toluene_xyz_file,
    ):
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "sterimol",
                "-f",
                sterimol_toluene_xyz_file,
                "--axis",
                "2",
                "3",
                "--auto-fragment",
            ],
        )
        assert result.exit_code == 1
        assert "is in a ring" in result.output
        assert not list(tmp_path.glob("chemsmart_*"))
