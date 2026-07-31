"""
Tests for the ``mol`` CLI subcommands (PyMOL visualization jobs).

Each test uses :class:`click.testing.CliRunner` to invoke the ``mol``
group and :mod:`unittest.mock` to intercept the job constructor so that
the merged settings can be inspected without running an actual PyMOL job.
"""

import os
import shutil
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.mol.mol import mol


def run_mol_and_capture_kwargs(job_class_path, cli_args):
    """Run the ``mol`` CLI with a patched job class and capture call kwargs."""
    runner = CliRunner()
    with patch(job_class_path) as mock_job_cls:
        mock_job_cls.return_value = MagicMock()
        result = runner.invoke(
            mol,
            cli_args,
            obj={},
            catch_exceptions=False,
        )
        call_kwargs = None
        if mock_job_cls.call_args is not None:
            call_kwargs = mock_job_cls.call_args.kwargs
    return result, call_kwargs


class TestMolCLINciCommand:
    """CLI tests for the ``nci`` subcommand."""

    def test_nci_job_creation(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.nci.PyMOLNCIJob",
            ["-f", single_molecule_xyz_file, "nci"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs is not None, "PyMOLNCIJob was never instantiated"

    def test_nci_coordinates_option(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.nci.PyMOLNCIJob",
            [
                "-f",
                single_molecule_xyz_file,
                "nci",
                "-c",
                "[[1,2],[3,4]]",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["coordinates"] == [[1, 2], [3, 4]]

    def test_nci_invalid_coordinates_raises(self, single_molecule_xyz_file):
        with pytest.raises(ValueError, match="Invalid coordinates input"):
            run_mol_and_capture_kwargs(
                "chemsmart.jobs.mol.nci.PyMOLNCIJob",
                [
                    "-f",
                    single_molecule_xyz_file,
                    "nci",
                    "-c",
                    "not-a-literal(",
                ],
            )


class TestMolCLIMoCommand:
    """CLI tests for the ``mo`` (molecular orbital) subcommand."""

    def test_mo_job_creation(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.mo.PyMOLMOJob",
            ["-f", single_molecule_xyz_file, "mo"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs is not None, "PyMOLMOJob was never instantiated"

    def test_mo_homo_flag(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.mo.PyMOLMOJob",
            ["-f", single_molecule_xyz_file, "mo", "--homo"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["homo"] is True


class TestMolCLIMovieCommand:
    """CLI tests for the ``movie`` subcommand."""

    def test_movie_job_creation(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.movie.PyMOLMovieJob",
            ["-f", single_molecule_xyz_file, "movie"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs is not None, "PyMOLMovieJob was never instantiated"

    def test_movie_overwrite_flag(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.movie.PyMOLMovieJob",
            ["-f", single_molecule_xyz_file, "movie", "-o"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["overwrite"] is True


class TestMolCLIIrcCommand:
    """CLI tests for the ``irc`` subcommand."""

    def test_irc_full_run_uses_all_file(self):
        runner = CliRunner()
        with patch("chemsmart.jobs.mol.irc.PyMOLIRCMovieJob") as mock_job_cls:
            mock_job_cls.from_files.return_value = MagicMock()
            result = runner.invoke(
                mol,
                ["irc", "-a", "full_irc.log"],
                obj={},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_job_cls.from_files.assert_called_once()
        assert (
            mock_job_cls.from_files.call_args.kwargs["all_file"]
            == "full_irc.log"
        )

    def test_irc_reactant_and_product_files(self):
        runner = CliRunner()
        with patch("chemsmart.jobs.mol.irc.PyMOLIRCMovieJob") as mock_job_cls:
            mock_job_cls.from_files.return_value = MagicMock()
            result = runner.invoke(
                mol,
                [
                    "irc",
                    "-r",
                    "reactant.log",
                    "-p",
                    "product.log",
                ],
                obj={},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert (
            mock_job_cls.from_files.call_args.kwargs["reactant_file"]
            == "reactant.log"
        )
        assert (
            mock_job_cls.from_files.call_args.kwargs["product_file"]
            == "product.log"
        )


class TestMolCLISpinCommand:
    """CLI tests for the ``spin`` subcommand."""

    def test_spin_job_creation(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.spin.PyMOLSpinJob",
            ["-f", single_molecule_xyz_file, "spin"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs is not None, "PyMOLSpinJob was never instantiated"

    def test_spin_basename_uses_provided_label(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.spin.PyMOLSpinJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-l",
                "custom_label",
                "spin",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["spin_basename"] == "custom_label"
        assert kwargs["label"] == "custom_label"


class TestMolCLIAlignCommand:
    """CLI tests for the ``align`` subcommand."""

    def test_align_requires_input_files(self):
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["align"],
        )
        assert result.exit_code != 0

    def test_align_multi_structure_single_file(
        self, multiple_molecules_xyz_file
    ):
        """A single multi-structure file with default index (``:``) aligns all structures."""
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["-f", multiple_molecules_xyz_file, "align"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs is not None, "PyMOLAlignJob was never instantiated"
        assert len(kwargs["molecule"]) >= 2

    def test_align_directory_without_filetype_raises(
        self, structure_test_directory
    ):
        """Directory given via the -p/--program group path (so
        ctx.obj["directory"] is set but ctx.obj["filetype"] is not)
        reaches align's own "directory without filetype" guard."""
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            [
                "-d",
                os.path.join(structure_test_directory, "xyz"),
                "-p",
                "gaussian",
                "align",
            ],
        )
        assert result.exit_code != 0
        assert "no filetype provided" in result.output

    def test_align_directory_with_filetype(
        self,
        tmp_path,
        single_molecule_xyz_file,
        two_rotated_molecules_xyz_file,
    ):
        """Directory + filetype: align loads and aligns one structure
        from each matched file in the directory."""
        shutil.copy(single_molecule_xyz_file, tmp_path / "a.xyz")
        shutil.copy(two_rotated_molecules_xyz_file, tmp_path / "b.xyz")
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["-d", str(tmp_path), "-t", "xyz", "align"],
        )
        assert result.exit_code == 0, result.output
        assert len(kwargs["molecule"]) == 2

    def test_align_multiple_files_two_structures_label(
        self, single_molecule_xyz_file, two_rotated_molecules_xyz_file
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-f",
                two_rotated_molecules_xyz_file,
                "align",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(kwargs["molecule"]) == 2
        assert kwargs["label"].endswith("_and_1_structure_align")

    def test_align_multiple_files_more_than_two_structures_label(
        self,
        single_molecule_xyz_file,
        two_rotated_molecules_xyz_file,
        visualized_1_mer_xyz_file,
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-f",
                two_rotated_molecules_xyz_file,
                "-f",
                visualized_1_mer_xyz_file,
                "align",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(kwargs["molecule"]) == 3
        assert kwargs["label"].endswith("_and_2_structures_align")

    def test_align_explicit_label_used_directly(
        self, single_molecule_xyz_file, two_rotated_molecules_xyz_file
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-f",
                two_rotated_molecules_xyz_file,
                "-l",
                "mycustom",
                "align",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"] == "mycustom"

    def test_align_not_enough_molecules_raises(self, single_molecule_xyz_file):
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["-f", single_molecule_xyz_file, "-i", "1", "align"],
        )
        assert result.exit_code != 0
        assert "at least 2 molecules" in result.output
        assert "may not select enough structures" in result.output

    def test_align_directory_no_files_found_for_filetype(
        self, tmp_path, single_molecule_xyz_file
    ):
        shutil.copy(single_molecule_xyz_file, tmp_path / "a.xyz")
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["-d", str(tmp_path), "-t", "nonexistentext", "align"],
        )
        assert result.exit_code != 0
        assert "No files found with extension" in result.output

    def test_align_directory_out_of_range_index_raises_bad_parameter(
        self,
        tmp_path,
        single_molecule_xyz_file,
        two_rotated_molecules_xyz_file,
    ):
        """Covers the per-file ValueError -> click.BadParameter wrapping
        shared by both the directory and multi-filenames branches."""
        shutil.copy(single_molecule_xyz_file, tmp_path / "a.xyz")
        shutil.copy(two_rotated_molecules_xyz_file, tmp_path / "b.xyz")
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["-d", str(tmp_path), "-t", "xyz", "-i", "5", "align"],
        )
        assert result.exit_code != 0
        assert "Error processing file" in result.output
        assert "out of range" in result.output

    def test_align_single_file_out_of_range_index_raises_bad_parameter(
        self, single_molecule_xyz_file
    ):
        """Covers the single-file branch's own ValueError ->
        click.BadParameter wrapping (distinct from the per-file helper
        used by the directory/multi-filenames branches)."""
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            ["-f", single_molecule_xyz_file, "-i", "100", "align"],
        )
        assert result.exit_code != 0
        assert "out of range" in result.output
        assert "Error processing file" not in result.output

    def test_align_multiple_files_out_of_range_index_raises_bad_parameter(
        self, single_molecule_xyz_file, two_rotated_molecules_xyz_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.align.PyMOLAlignJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-f",
                two_rotated_molecules_xyz_file,
                "-i",
                "5",
                "align",
            ],
        )
        assert result.exit_code != 0
        assert "Error processing file" in result.output
        assert "out of range" in result.output


class TestMolCLIVisualizeCommand:
    """CLI tests for the ``visualize`` subcommand."""

    def test_visualize_basic_job_creation(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLVisualizationJob",
            ["-f", single_molecule_xyz_file, "visualize"],
        )
        assert result.exit_code == 0, result.output
        assert (
            kwargs is not None
        ), "PyMOLVisualizationJob was never instantiated"

    def test_visualize_hybrid_job_creation(self, single_molecule_xyz_file):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLHybridVisualizationJob",
            [
                "-f",
                single_molecule_xyz_file,
                "visualize",
                "-H",
                "-G",
                "1,2-5",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (
            kwargs is not None
        ), "PyMOLHybridVisualizationJob was never instantiated"
        assert kwargs["groups"] == ("1,2-5",)

    def test_visualize_derived_style_job_creation(
        self, single_molecule_xyz_file
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLScientificStyleVisualizationJob",
            [
                "-f",
                single_molecule_xyz_file,
                "visualize",
                "-s",
                "editorial-minimal",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (
            kwargs is not None
        ), "PyMOLScientificStyleVisualizationJob was never instantiated"
        assert kwargs["style"] == "editorial_minimal"

    def test_visualize_hybrid_with_derived_style_raises(
        self, single_molecule_xyz_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLHybridVisualizationJob",
            [
                "-f",
                single_molecule_xyz_file,
                "visualize",
                "-H",
                "-s",
                "editorial-minimal",
            ],
        )
        assert result.exit_code != 0
