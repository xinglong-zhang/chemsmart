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


NCI_JOB = "chemsmart.jobs.mol.nci.PyMOLNCIJob"


class TestMolCLIGroupValidation:
    """Direct CLI tests for the ``mol`` group callback's own parsing
    and validation logic (chemsmart:cli/mol/mol.py), using the ``nci``
    subcommand as a lightweight vehicle since group-level behavior is
    independent of which subcommand is ultimately invoked."""

    def test_index_and_structure_index_mutually_exclusive(
        self, single_molecule_xyz_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            NCI_JOB,
            [
                "-f",
                single_molecule_xyz_file,
                "-i",
                "1",
                "--si",
                "1",
                "nci",
            ],
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_structure_index_used_as_index(self, multiple_molecules_xyz_file):
        """``--si`` alone is treated as an alias for ``-i``."""
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB,
            ["-f", multiple_molecules_xyz_file, "--si", "1", "nci"],
        )
        assert result.exit_code == 0, result.output
        assert len(kwargs["molecule"]) == 1

    def test_chemsmart_db_requires_exactly_one_selector(
        self, database_chemsmart_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            NCI_JOB, ["-f", database_chemsmart_file, "nci"]
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_rejects_multiple_selectors(
        self, database_chemsmart_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            NCI_JOB,
            [
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "--sid",
                "f751bb2c27e2",
                "nci",
            ],
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_index_requires_record_selector(
        self, database_chemsmart_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            NCI_JOB,
            [
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "-i",
                "1",
                "nci",
            ],
        )
        assert result.exit_code != 0
        assert "can only be used together with" in result.output

    def test_directory_and_filetype_uses_auto_label(self, tmp_path):
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB, ["-d", str(tmp_path), "-t", "xyz", "nci"]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"] == f"all_xyz_files_in_{tmp_path.name}"

    def test_directory_and_filetype_respects_explicit_label(self, tmp_path):
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB, ["-d", str(tmp_path), "-t", "xyz", "-l", "custom", "nci"]
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"] == "custom"

    def test_directory_and_program_uses_auto_label(self, tmp_path):
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB, ["-d", str(tmp_path), "-p", "gaussian", "nci"]
        )
        assert result.exit_code == 0, result.output
        assert (
            kwargs["label"]
            == f"all_output_files_from_gaussian_in_{tmp_path.name}"
        )

    def test_directory_and_program_respects_explicit_label(self, tmp_path):
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB,
            ["-d", str(tmp_path), "-p", "gaussian", "-l", "custom", "nci"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"] == "custom"

    def test_filenames_and_pubchem_both_specified_raises(
        self, single_molecule_xyz_file
    ):
        with pytest.raises(ValueError, match="have been specified"):
            run_mol_and_capture_kwargs(
                NCI_JOB,
                ["-f", single_molecule_xyz_file, "-P", "water", "nci"],
            )

    def test_multiple_filenames_non_align_task_raises(
        self, single_molecule_xyz_file, two_rotated_molecules_xyz_file
    ):
        with pytest.raises(ValueError, match="can only process one file"):
            run_mol_and_capture_kwargs(
                NCI_JOB,
                [
                    "-f",
                    single_molecule_xyz_file,
                    "-f",
                    two_rotated_molecules_xyz_file,
                    "nci",
                ],
            )

    def test_label_and_append_label_both_given_raises(
        self, single_molecule_xyz_file
    ):
        with pytest.raises(ValueError, match="not both"):
            run_mol_and_capture_kwargs(
                NCI_JOB,
                [
                    "-f",
                    single_molecule_xyz_file,
                    "-l",
                    "x",
                    "-a",
                    "y",
                    "nci",
                ],
            )

    def test_pubchem_only_without_label_crashes(self):
        """Regression test: -P/--pubchem given without -l/--label
        crashes, since the default-label derivation unconditionally
        does os.path.basename(filenames) even though filenames stays
        None for a pubchem-only invocation. Mirrors the same
        pubchem-only-without-label bug already documented for the
        gaussian/orca CLIs (BUGS_FOUND.md #25)."""
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[MagicMock()],
        ):
            with pytest.raises(TypeError):
                run_mol_and_capture_kwargs(NCI_JOB, ["-P", "water", "nci"])

    def test_append_label_non_chemsmart_db_file(
        self, single_molecule_xyz_file
    ):
        """The is_chemsmart_db-specific suffix branches (SID/RID/RI/MID)
        are skipped entirely for an ordinary (non-database) file."""
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB,
            ["-f", single_molecule_xyz_file, "-a", "suffix", "nci"],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"] == "crest_best_suffix"

    def test_pubchem_with_label_succeeds(self):
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[MagicMock()],
        ) as mock_from_pubchem:
            result, kwargs = run_mol_and_capture_kwargs(
                NCI_JOB, ["-P", "water", "-l", "water_mol", "nci"]
            )
        assert result.exit_code == 0, result.output
        mock_from_pubchem.assert_called_once_with(
            identifier="water", return_list=True
        )
        assert kwargs["label"] == "water_mol"

    @pytest.mark.parametrize(
        "selector_flag,selector_value",
        [
            ("--ri", "1"),
            ("--rid", "6de213a0"),
            ("--sid", "f751bb2c27e2"),
            ("--mid", "YGQDVTLERYVWNN-UHFFFAOYSA-N"),
        ],
    )
    def test_chemsmart_db_default_label_selector_branches(
        self, database_chemsmart_file, selector_flag, selector_value
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB,
            [
                "-f",
                database_chemsmart_file,
                selector_flag,
                selector_value,
                "nci",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"].startswith("chemsmart_")

    @pytest.mark.parametrize(
        "selector_flag,selector_value",
        [
            ("--ri", "1"),
            ("--rid", "6de213a0"),
            ("--sid", "f751bb2c27e2"),
            ("--mid", "YGQDVTLERYVWNN-UHFFFAOYSA-N"),
        ],
    )
    def test_chemsmart_db_append_label_selector_branches(
        self, database_chemsmart_file, selector_flag, selector_value
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            NCI_JOB,
            [
                "-f",
                database_chemsmart_file,
                selector_flag,
                selector_value,
                "-a",
                "suffix",
                "nci",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["label"].endswith("_suffix")


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

    def test_mo_invalid_coordinates_raises(self, single_molecule_xyz_file):
        with pytest.raises(ValueError, match="Invalid coordinates input"):
            run_mol_and_capture_kwargs(
                "chemsmart.jobs.mol.mo.PyMOLMOJob",
                [
                    "-f",
                    single_molecule_xyz_file,
                    "mo",
                    "-c",
                    "not-a-literal(",
                ],
            )


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

    def test_movie_invalid_coordinates_raises(self, single_molecule_xyz_file):
        with pytest.raises(ValueError, match="Invalid coordinates input"):
            run_mol_and_capture_kwargs(
                "chemsmart.jobs.mol.movie.PyMOLMovieJob",
                [
                    "-f",
                    single_molecule_xyz_file,
                    "movie",
                    "-c",
                    "not-a-literal(",
                ],
            )


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

    def test_irc_invalid_coordinates_raises(self):
        runner = CliRunner()
        with patch("chemsmart.jobs.mol.irc.PyMOLIRCMovieJob") as mock_job_cls:
            mock_job_cls.from_files.return_value = MagicMock()
            with pytest.raises(ValueError, match="Invalid coordinates input"):
                runner.invoke(
                    mol,
                    ["irc", "-a", "full_irc.log", "-c", "not-a-literal("],
                    obj={},
                    catch_exceptions=False,
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

    def test_spin_invalid_coordinates_raises(self, single_molecule_xyz_file):
        with pytest.raises(ValueError, match="Invalid coordinates input"):
            run_mol_and_capture_kwargs(
                "chemsmart.jobs.mol.spin.PyMOLSpinJob",
                [
                    "-f",
                    single_molecule_xyz_file,
                    "spin",
                    "-c",
                    "not-a-literal(",
                ],
            )


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

    def test_visualize_invalid_coordinates_raises(
        self, single_molecule_xyz_file
    ):
        with pytest.raises(ValueError, match="Invalid coordinates input"):
            run_mol_and_capture_kwargs(
                "chemsmart.jobs.mol.visualize.PyMOLVisualizationJob",
                [
                    "-f",
                    single_molecule_xyz_file,
                    "visualize",
                    "-c",
                    "not-a-literal(",
                ],
            )

    def test_visualize_hybrid_only_option_without_hybrid_raises(
        self, single_molecule_xyz_file
    ):
        result, _ = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLVisualizationJob",
            [
                "-f",
                single_molecule_xyz_file,
                "visualize",
                "--surface-color",
                "blue",
            ],
        )
        assert result.exit_code != 0
        assert "-H/--hybrid" in result.output

    def test_visualize_hybrid_only_option_with_hybrid_is_forwarded(
        self, single_molecule_xyz_file
    ):
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLHybridVisualizationJob",
            [
                "-f",
                single_molecule_xyz_file,
                "visualize",
                "-H",
                "--surface-color",
                "blue",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["surface_color"] == "blue"

    def test_visualize_non_derived_style_job_creation(
        self, single_molecule_xyz_file
    ):
        """A -s value that isn't one of the derived
        zhang_group_scientific_styles.py styles (e.g. "cylview") is
        passed straight through as the plain style kwarg."""
        result, kwargs = run_mol_and_capture_kwargs(
            "chemsmart.jobs.mol.visualize.PyMOLVisualizationJob",
            [
                "-f",
                single_molecule_xyz_file,
                "visualize",
                "-s",
                "cylview",
            ],
        )
        assert result.exit_code == 0, result.output
        assert kwargs["style"] == "cylview"


def _invoke_mol_qmmm_callback(**overrides):
    """Direct-invocation helper for BUGS_FOUND.md #38: `mol_qmmm` has
    no subcommands registered and is never attached to any parent CLI
    group, so it cannot be reached via CliRunner (Click raises "Missing
    command" before the group callback ever runs). This bypasses
    Click's command routing to unit-test the callback body itself."""
    import click

    from chemsmart.cli.mol.mol import mol_qmmm

    kwargs = dict(
        filenames=None,
        label=None,
        append_label=None,
        index=None,
        directory=None,
        filetype=None,
        program=None,
        pubchem=None,
    )
    kwargs.update(overrides)
    ctx = click.Context(mol_qmmm, obj={})
    ctx.invoke(mol_qmmm.callback, **kwargs)
    return ctx.obj


class TestMolQmmmGroupDirectInvocation:
    """See BUGS_FOUND.md #38: mol_qmmm is orphaned code with no
    reachable entry point. These tests exercise its callback body
    directly via ctx.invoke to get unit coverage on logic that no real
    CLI invocation can ever reach."""

    def test_directory_and_filetype(self, tmp_path, single_molecule_xyz_file):
        shutil.copy(single_molecule_xyz_file, tmp_path / "a.xyz")
        obj = _invoke_mol_qmmm_callback(
            directory=str(tmp_path), filetype="xyz"
        )
        assert obj["label"] == f"all_xyz_files_in_{tmp_path.name}"
        assert len(obj["molecules"]) == 1

    def test_directory_and_filetype_respects_explicit_label(
        self, tmp_path, single_molecule_xyz_file
    ):
        shutil.copy(single_molecule_xyz_file, tmp_path / "a.xyz")
        obj = _invoke_mol_qmmm_callback(
            directory=str(tmp_path), filetype="xyz", label="custom"
        )
        assert obj["label"] == "custom"

    def test_directory_and_program(self, tmp_path):
        """An empty directory yields zero matched output files, so no
        real Gaussian/ORCA log parsing is needed to exercise this
        branch (mirrors the equivalent `mol`-group test)."""
        obj = _invoke_mol_qmmm_callback(
            directory=str(tmp_path), program="gaussian"
        )
        assert (
            obj["label"]
            == f"all_output_files_from_gaussian_in_{tmp_path.name}"
        )
        assert obj["molecules"] == []

    def test_no_filenames_no_pubchem_warns_and_returns_none(self):
        obj = _invoke_mol_qmmm_callback()
        assert obj["molecules"] is None
        assert obj["label"] is None
        assert obj["qmmm"] is True

    def test_filenames_and_pubchem_both_specified_raises(self):
        with pytest.raises(ValueError, match="have been specified"):
            _invoke_mol_qmmm_callback(filenames=("a.xyz",), pubchem="water")

    def test_single_filename_loads_qmmm_molecule(
        self, single_molecule_xyz_file
    ):
        obj = _invoke_mol_qmmm_callback(filenames=(single_molecule_xyz_file,))
        assert obj["label"] == "crest_best"
        assert len(obj["molecules"]) == 1

    def test_multiple_filenames_early_return(self):
        obj = _invoke_mol_qmmm_callback(filenames=("a.xyz", "b.xyz"))
        assert obj["molecules"] is None
        assert obj["filenames"] == ("a.xyz", "b.xyz")

    def test_pubchem_success(self):
        with patch(
            "chemsmart.io.molecules.structure.QMMMMolecule.from_pubchem",
            return_value=[MagicMock()],
        ) as mock_from_pubchem:
            obj = _invoke_mol_qmmm_callback(pubchem="water", label="water_mol")
        assert obj["label"] == "water_mol"
        mock_from_pubchem.assert_called_once_with(
            identifier="water", return_list=True
        )

    def test_label_and_append_label_both_given_raises(
        self, single_molecule_xyz_file
    ):
        with pytest.raises(ValueError, match="not both"):
            _invoke_mol_qmmm_callback(
                filenames=(single_molecule_xyz_file,),
                label="x",
                append_label="y",
            )

    def test_append_label_suffix(self, single_molecule_xyz_file):
        obj = _invoke_mol_qmmm_callback(
            filenames=(single_molecule_xyz_file,), append_label="suffix"
        )
        assert obj["label"] == "crest_best_suffix"

    def test_index_selects_specific_structure(
        self, multiple_molecules_xyz_file
    ):
        obj = _invoke_mol_qmmm_callback(
            filenames=(multiple_molecules_xyz_file,), index="2"
        )
        assert len(obj["molecules"]) == 1

    def test_default_index_uses_last_molecule(
        self, multiple_molecules_xyz_file
    ):
        obj = _invoke_mol_qmmm_callback(
            filenames=(multiple_molecules_xyz_file,)
        )
        assert len(obj["molecules"]) == 1

    def test_molecules_converted_to_qmmm_molecule(
        self, single_molecule_xyz_file
    ):
        from chemsmart.io.molecules.structure import QMMMMolecule

        obj = _invoke_mol_qmmm_callback(filenames=(single_molecule_xyz_file,))
        assert all(isinstance(m, QMMMMolecule) for m in obj["molecules"])
