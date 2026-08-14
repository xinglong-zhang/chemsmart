"""Tests for chemsmart.cli.grouper (the grouper CLI group and subcommands)."""

import importlib
from unittest.mock import MagicMock

import click
import pytest
from click.testing import CliRunner

from chemsmart.cli.grouper import grouper
from chemsmart.cli.grouper.grouper import (
    THERMOCHEMISTRY_DEFAULT_KWARGS,
    _extract_conformer_id,
    _extract_energy_based_on_energy_type,
    _extract_last_energy_from_output_file,
    _find_matching_sp_file,
    _get_label,
    _load_molecules_from_directory,
    build_thermochemistry_kwargs,
    create_grouper_job_from_context,
)

# The `chemsmart.cli.grouper` package's __init__.py does
# ``from .grouper import grouper``, which rebinds the package attribute
# ``grouper`` to the Click group object rather than the submodule.
# String-based ``mocker.patch("chemsmart.cli.grouper.grouper.X")`` resolves
# via attribute traversal and would therefore hit the Click group, not the
# submodule. Look the submodule up via ``sys.modules`` (through
# ``importlib.import_module``, which checks the module cache directly)
# and patch attributes on that object instead.
grouper_module = importlib.import_module("chemsmart.cli.grouper.grouper")


def invoke_grouper(args, mocker, mock_job_return=None):
    mock_job_cls = mocker.patch("chemsmart.jobs.grouper.GrouperJob")
    if mock_job_return is not None:
        mock_job_cls.return_value = mock_job_return
    runner = CliRunner()
    result = runner.invoke(grouper, args, catch_exceptions=False)
    return result, mock_job_cls


class TestBuildThermochemistryKwargs:
    def test_defaults_when_nothing_given(self):
        kwargs = build_thermochemistry_kwargs()
        for key, value in THERMOCHEMISTRY_DEFAULT_KWARGS.items():
            assert kwargs[key] == value

    def test_both_cutoffs_raises(self):
        with pytest.raises(ValueError, match="Cannot specify both"):
            build_thermochemistry_kwargs(
                cutoff_entropy_grimme=100, cutoff_entropy_truhlar=110
            )

    def test_grimme_cutoff_sets_entropy_method(self):
        kwargs = build_thermochemistry_kwargs(cutoff_entropy_grimme=100)
        assert kwargs["s_freq_cutoff"] == 100
        assert kwargs["entropy_method"] == "grimme"

    def test_truhlar_cutoff_sets_entropy_method(self):
        kwargs = build_thermochemistry_kwargs(cutoff_entropy_truhlar=110)
        assert kwargs["s_freq_cutoff"] == 110
        assert kwargs["entropy_method"] == "truhlar"

    def test_overrides_applied(self):
        kwargs = build_thermochemistry_kwargs(
            cutoff_enthalpy=50,
            concentration=2.0,
            pressure=0.5,
            temperature=300.0,
            alpha=8,
            weighted=False,
            energy_units="eV",
            check_imaginary_frequencies=False,
        )
        assert kwargs["h_freq_cutoff"] == 50
        assert kwargs["concentration"] == 2.0
        assert kwargs["pressure"] == 0.5
        assert kwargs["temperature"] == 300.0
        assert kwargs["alpha"] == 8
        assert kwargs["use_weighted_mass"] is False
        assert kwargs["energy_units"] == "eV"
        assert kwargs["check_imaginary_frequencies"] is False


class TestExtractConformerId:
    def test_matches_underscore_pattern(self):
        assert _extract_conformer_id("structure1_c12_opt.log") == "c12"

    def test_matches_dot_pattern(self):
        assert _extract_conformer_id("structure1_c5.log") == "c5"

    def test_no_match_returns_none(self):
        assert _extract_conformer_id("structure1_opt.log") is None


class TestGetLabel:
    def test_both_label_and_append_label_raises(self):
        with pytest.raises(ValueError, match="but not both"):
            _get_label("mylabel", "suffix", "base")

    def test_append_label_appends_to_base(self):
        assert _get_label(None, "suffix", "base") == "base_suffix"

    def test_label_returns_label(self):
        assert _get_label("mylabel", None, "base") == "mylabel"

    def test_neither_returns_base(self):
        assert _get_label(None, None, "base") == "base"


class TestExtractEnergyBasedOnEnergyType:
    def _mock_thermo(self, **output_attrs):
        thermo = MagicMock()
        output = MagicMock()
        for key, value in output_attrs.items():
            setattr(output, key, value)
        thermo.file_object = output
        return thermo

    def test_energy_type_e_returns_last_energy(self):
        thermo = self._mock_thermo(energies=[1.0, 2.0, 3.0])
        assert _extract_energy_based_on_energy_type(thermo, "E") == 3.0

    def test_energy_type_e_no_energies_returns_none(self):
        thermo = self._mock_thermo(energies=[])
        assert _extract_energy_based_on_energy_type(thermo, "E") is None

    def test_energy_type_h_returns_enthalpy(self):
        thermo = self._mock_thermo(enthalpy=-5.0)
        assert _extract_energy_based_on_energy_type(thermo, "H") == -5.0

    def test_energy_type_h_attribute_error_returns_none(self):
        output = MagicMock(spec=[])
        thermo = MagicMock()
        thermo.file_object = output
        assert _extract_energy_based_on_energy_type(thermo, "H") is None

    def test_energy_type_g_returns_gibbs(self):
        thermo = self._mock_thermo(gibbs_free_energy=-6.0)
        assert _extract_energy_based_on_energy_type(thermo, "G") == -6.0

    def test_energy_type_g_attribute_error_returns_none(self):
        output = MagicMock(spec=[])
        thermo = MagicMock()
        thermo.file_object = output
        assert _extract_energy_based_on_energy_type(thermo, "G") is None

    def test_energy_type_qhh_uses_qrrho_enthalpy(self, mocker):
        mocker.patch.object(
            grouper_module,
            "energy_conversion",
            return_value=-7.0,
        )
        thermo = MagicMock()
        thermo.qrrho_enthalpy = -12345.0
        result = _extract_energy_based_on_energy_type(thermo, "QHH")
        assert result == -7.0

    def test_energy_type_qhh_falls_back_to_h_on_attribute_error(self, mocker):
        thermo = MagicMock(spec=["file_object"])
        thermo.file_object = MagicMock(enthalpy=-5.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHH")
        assert result == -5.0

    def test_energy_type_qhh_falls_back_to_h_when_none(self):
        thermo = MagicMock()
        thermo.qrrho_enthalpy = None
        thermo.file_object = MagicMock(enthalpy=-5.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHH")
        assert result == -5.0

    def test_energy_type_qhg_uses_qrrho_gibbs(self, mocker):
        mocker.patch.object(
            grouper_module,
            "energy_conversion",
            return_value=-8.0,
        )
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = -54321.0
        result = _extract_energy_based_on_energy_type(thermo, "QHG")
        assert result == -8.0

    def test_energy_type_qhg_falls_back_to_g_on_attribute_error(self):
        thermo = MagicMock(spec=["file_object"])
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHG")
        assert result == -6.0

    def test_energy_type_qhg_falls_back_to_g_when_none(self):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0)
        result = _extract_energy_based_on_energy_type(thermo, "QHG")
        assert result == -6.0

    def test_sp_qhg_full_success(self, mocker):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        thermo.filename = "/tmp/foo_c1_opt.log"
        mocker.patch.object(
            grouper_module,
            "_find_matching_sp_file",
            return_value="/tmp/sp/foo_c1_sp.log",
        )
        mocker.patch.object(
            grouper_module,
            "_extract_last_energy_from_output_file",
            return_value=-4.0,
        )
        result = _extract_energy_based_on_energy_type(thermo, "SP_QHG")
        assert result == -6.0 - (-5.0) + (-4.0)

    def test_sp_qhg_no_filename_returns_qh_gibbs(self):
        thermo = MagicMock(spec=["file_object", "qrrho_gibbs_free_energy"])
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        result = _extract_energy_based_on_energy_type(thermo, "SP_QHG")
        assert result == -6.0

    def test_sp_qhg_no_matching_sp_file_raises(self, mocker):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        thermo.filename = "/tmp/foo_c1_opt.log"
        mocker.patch.object(
            grouper_module,
            "_find_matching_sp_file",
            return_value=None,
        )
        with pytest.raises(click.ClickException, match="matching SP file"):
            _extract_energy_based_on_energy_type(thermo, "SP_QHG")

    def test_sp_qhg_esolv_none_raises(self, mocker):
        thermo = MagicMock()
        thermo.qrrho_gibbs_free_energy = None
        thermo.file_object = MagicMock(gibbs_free_energy=-6.0, energies=[-5.0])
        thermo.filename = "/tmp/foo_c1_opt.log"
        mocker.patch.object(
            grouper_module,
            "_find_matching_sp_file",
            return_value="/tmp/sp/foo_c1_sp.log",
        )
        mocker.patch.object(
            grouper_module,
            "_extract_last_energy_from_output_file",
            return_value=None,
        )
        with pytest.raises(click.ClickException, match="sp_qh_gibbs"):
            _extract_energy_based_on_energy_type(thermo, "SP_QHG")

    def test_unsupported_energy_type_raises(self):
        thermo = MagicMock()
        with pytest.raises(click.ClickException, match="Unsupported"):
            _extract_energy_based_on_energy_type(thermo, "BOGUS")


class TestFindMatchingSpFile:
    def test_raises_when_sp_folder_missing(self, tmp_path):
        source = tmp_path / "foo_c1_opt.log"
        source.write_text("")
        with pytest.raises(ValueError, match="Can't find SP output folder"):
            _find_matching_sp_file(str(source))

    def test_finds_matching_suffix_file(self, tmp_path):
        sp_folder = tmp_path / "sp"
        sp_folder.mkdir()
        (sp_folder / "foo_c1_opt_sp.log").write_text("")
        source = tmp_path / "foo_c1_opt.log"
        source.write_text("")
        result = _find_matching_sp_file(str(source))
        assert result == str(sp_folder / "foo_c1_opt_sp.log")

    def test_returns_none_when_no_candidates(self, tmp_path):
        sp_folder = tmp_path / "sp"
        sp_folder.mkdir()
        source = tmp_path / "foo_c1_opt.log"
        source.write_text("")
        assert _find_matching_sp_file(str(source)) is None


class TestExtractLastEnergyFromOutputFile:
    def test_gaussian_program_returns_last_energy(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="gaussian",
        )
        mock_output = MagicMock(normal_termination=True, energies=[-1.0, -2.0])
        mocker.patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        )
        result = _extract_last_energy_from_output_file("foo.log")
        assert result == -2.0

    def test_orca_program_returns_last_energy(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="orca",
        )
        mock_output = MagicMock(normal_termination=True, energies=[-3.0, -4.0])
        mocker.patch(
            "chemsmart.io.orca.output.ORCAOutput",
            return_value=mock_output,
        )
        result = _extract_last_energy_from_output_file("foo.out")
        assert result == -4.0

    def test_unknown_program_returns_none(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value=None,
        )
        assert _extract_last_energy_from_output_file("foo.xyz") is None

    def test_not_normal_termination_returns_none(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="gaussian",
        )
        mock_output = MagicMock(normal_termination=False)
        mocker.patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        )
        assert _extract_last_energy_from_output_file("foo.log") is None

    def test_no_energies_returns_none(self, mocker):
        mocker.patch.object(
            grouper_module,
            "get_program_type_from_file",
            return_value="gaussian",
        )
        mock_output = MagicMock(normal_termination=True, energies=[])
        mocker.patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        )
        assert _extract_last_energy_from_output_file("foo.log") is None


class TestGrouperCliValidation:
    def test_both_filenames_and_directory_raises(
        self, mocker, multiple_molecules_xyz_file, tmp_path
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-d",
                str(tmp_path),
                "-p",
                "gaussian",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "Cannot specify both" in result.output

    def test_directory_without_program_or_filetype_raises(
        self, mocker, tmp_path
    ):
        result, _ = invoke_grouper(
            ["-d", str(tmp_path), "rmsd"],
            mocker,
        )
        assert result.exit_code != 0
        assert "Must specify -p/--program" in result.output

    def test_multiple_filenames_rejected(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-f",
                multiple_molecules_xyz_file,
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "only single file input" in result.output

    def test_non_e_energy_type_for_xyz_rejected(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-E",
                "H",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "Only energy" in result.output

    def test_no_input_specified_warns_and_returns(self, mocker):
        result, mock_job_cls = invoke_grouper(["rmsd"], mocker)
        assert result.exit_code == 0
        mock_job_cls.assert_called_once()
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["molecules"] is None


class TestGrouperCliSuccessPaths:
    def test_single_file_rmsd_success(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, "rmsd"], mocker
        )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["grouping_strategy"] == "rmsd"
        assert len(call_kwargs["molecules"]) > 1

    def test_label_option_used_as_grouper_label(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-l",
                "mylabel",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["label"] == "mylabel_rmsd"

    def test_directory_mode_success(self, mocker, tmp_path):
        mock_molecule = MagicMock(name="a")
        mocker.patch.object(
            grouper_module,
            "_load_molecules_from_directory",
            return_value=(
                [mock_molecule, mock_molecule],
                ["c1", "c2"],
                [],
                "auto_label",
            ),
        )
        result, mock_job_cls = invoke_grouper(
            ["-d", str(tmp_path), "-p", "gaussian", "rmsd"], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["conformer_ids"] == ["c1", "c2"]
        assert call_kwargs["label"] == "auto_label_rmsd"

    def test_fewer_than_two_molecules_warns_but_continues(
        self, mocker, tmp_path
    ):
        single_mol_file = tmp_path / "single.xyz"
        single_mol_file.write_text("1\ncomment\nH 0.0 0.0 0.0\n")
        result, mock_job_cls = invoke_grouper(
            ["-f", str(single_mol_file), "rmsd"], mocker
        )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()

    def test_tanimoto_default_fingerprint(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, "tanimoto"], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["fingerprint_type"] == "rdkit"

    def test_tanimoto_custom_fingerprint(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "tanimoto",
                "--fingerprint-type",
                "morgan",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["fingerprint_type"] == "morgan"

    def test_tfd_options(self, mocker, multiple_molecules_xyz_file):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "tfd",
                "--no-use-weights",
                "--max-dev",
                "spec",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["use_weights"] is False
        assert call_kwargs["max_dev"] == "spec"

    def test_irmsd_inversion_option(self, mocker, multiple_molecules_xyz_file):
        result, mock_job_cls = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "irmsd",
                "--inversion",
                "on",
            ],
            mocker,
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["inversion"] == "on"

    @pytest.mark.parametrize(
        "subcommand", ["hrmsd", "spyrmsd", "pymolrmsd", "energy"]
    )
    def test_simple_subcommands_succeed(
        self, mocker, multiple_molecules_xyz_file, subcommand
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, subcommand], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["grouping_strategy"] == subcommand


class TestGrouperCliStructuralSubcommands:
    @pytest.mark.parametrize(
        "subcommand,strategy",
        [
            ("connectivity", "connectivity"),
            ("formula", "formula"),
            ("isomorphism", "isomorphism"),
        ],
    )
    def test_succeeds_without_threshold_or_num_groups(
        self, mocker, multiple_molecules_xyz_file, subcommand, strategy
    ):
        result, mock_job_cls = invoke_grouper(
            ["-f", multiple_molecules_xyz_file, subcommand], mocker
        )
        assert result.exit_code == 0, result.output
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["grouping_strategy"] == strategy

    @pytest.mark.parametrize(
        "subcommand", ["connectivity", "formula", "isomorphism"]
    )
    def test_rejects_threshold_option(
        self, mocker, multiple_molecules_xyz_file, subcommand
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-T",
                "0.5",
                subcommand,
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "does not support threshold" in result.output

    @pytest.mark.parametrize(
        "subcommand", ["connectivity", "formula", "isomorphism"]
    )
    def test_rejects_num_groups_option(
        self, mocker, multiple_molecules_xyz_file, subcommand
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-N",
                "3",
                subcommand,
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "does not support num_groups" in result.output


class TestGrouperCliBothCutoffsRejected:
    def test_both_cutoffs_becomes_bad_parameter(
        self, mocker, multiple_molecules_xyz_file
    ):
        result, _ = invoke_grouper(
            [
                "-f",
                multiple_molecules_xyz_file,
                "-csg",
                "100",
                "-cst",
                "110",
                "rmsd",
            ],
            mocker,
        )
        assert result.exit_code != 0
        assert "Cannot specify both" in result.output


class TestFindMatchingSpFileNoSuffix:
    def test_source_without_suffix_skips_suffix_glob(self, tmp_path):
        sp_folder = tmp_path / "sp"
        sp_folder.mkdir()
        (sp_folder / "foo_c1_sp.log").write_text("")
        source = tmp_path / "foo_c1"
        source.write_text("")
        result = _find_matching_sp_file(str(source))
        assert result == str(sp_folder / "foo_c1_sp.log")


class TestLoadMoleculesFromDirectory:
    def test_missing_program_and_filetype_raises(self, mocker, tmp_path):
        mock_folder_cls = mocker.patch("chemsmart.io.folder.BaseFolder")
        mock_folder_cls.return_value = MagicMock()
        with pytest.raises(click.BadParameter, match="Must specify"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program=None, filetype=None
            )

    def test_no_output_files_found_raises(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            []
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        with pytest.raises(click.BadParameter, match="No matching output"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )

    def test_program_and_filetype_uses_combined_lookup(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_files_in_current_folder_by_program_and_suffix.return_value = (
            []
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        with pytest.raises(click.BadParameter, match="No matching output"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype="log"
            )
        mock_folder.get_all_files_in_current_folder_by_program_and_suffix.assert_called_once()

    def test_filetype_only_uses_suffix_lookup(self, mocker, tmp_path):
        mock_folder = MagicMock()
        mock_folder.get_all_files_in_current_folder_by_suffix.return_value = []
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        with pytest.raises(click.BadParameter, match="No matching output"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program=None, filetype="log"
            )

    def test_successful_load_with_patterned_and_unpatterned_files(
        self, mocker, tmp_path
    ):
        files = [
            str(tmp_path / "mol_c2_opt.log"),
            str(tmp_path / "mol_c1_opt.log"),
            str(tmp_path / "unpatterned.log"),
        ]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )

        mock_mol = MagicMock()
        mock_thermo = MagicMock()
        mock_thermo.molecule = mock_mol
        mock_thermo.file_object.energies = [-1.0]
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            return_value=mock_thermo,
        )

        molecules, conformer_ids, skipped_ids, common_label = (
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )
        )
        assert len(molecules) == 3
        # c1 sorted before c2, unpatterned file goes last.
        assert conformer_ids == ["c1", "c2", "unpatterned"]
        assert skipped_ids == []
        assert common_label == "mol"

    def test_skips_molecule_when_energy_extraction_fails(
        self, mocker, tmp_path
    ):
        files = [str(tmp_path / "mol_c1_opt.log")]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )

        mock_thermo = MagicMock()
        mock_thermo.molecule = MagicMock()
        mock_thermo.file_object.energies = []
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            return_value=mock_thermo,
        )

        with pytest.raises(click.BadParameter, match="No valid molecules"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )

    def test_skips_when_molecule_is_none(self, mocker, tmp_path):
        files = [str(tmp_path / "mol_c1_opt.log")]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )

        mock_thermo = MagicMock()
        mock_thermo.molecule = None
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            return_value=mock_thermo,
        )

        with pytest.raises(click.BadParameter, match="No valid molecules"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )

    def test_skips_on_thermochemistry_construction_error(
        self, mocker, tmp_path
    ):
        files = [str(tmp_path / "mol_c1_opt.log")]
        mock_folder = MagicMock()
        mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
            files
        )
        mocker.patch(
            "chemsmart.io.folder.BaseFolder", return_value=mock_folder
        )
        mocker.patch(
            "chemsmart.analysis.thermochemistry.Thermochemistry",
            side_effect=RuntimeError("boom"),
        )

        with pytest.raises(click.BadParameter, match="No valid molecules"):
            _load_molecules_from_directory(
                directory=str(tmp_path), program="gaussian", filetype=None
            )


class TestCreateGrouperJobFromContextThermoParameters:
    def test_qhg_energy_type_builds_thermo_parameters_string(self, mocker):
        mock_job_cls = mocker.patch("chemsmart.jobs.grouper.GrouperJob")
        ctx = MagicMock()
        ctx.obj = {
            "molecules": [MagicMock(), MagicMock()],
            "ignore_hydrogens": False,
            "num_procs": 1,
            "grouper_label": "mylabel",
            "num_groups": None,
            "conformer_ids": None,
            "skipped_ids": None,
            "energy_type": "qhG",
            "thermo_kwargs": {"temperature": 298.15, "alpha": None},
            "threshold": None,
            "matrix_format": "xlsx",
        }
        create_grouper_job_from_context(ctx, strategy="rmsd")
        _, call_kwargs = mock_job_cls.call_args
        assert "temperature=298.15" in call_kwargs["thermo_parameters"]

    def test_e_energy_type_leaves_thermo_parameters_none(self, mocker):
        mock_job_cls = mocker.patch("chemsmart.jobs.grouper.GrouperJob")
        ctx = MagicMock()
        ctx.obj = {
            "molecules": [MagicMock(), MagicMock()],
            "ignore_hydrogens": False,
            "num_procs": 1,
            "grouper_label": "mylabel",
            "num_groups": None,
            "conformer_ids": None,
            "skipped_ids": None,
            "energy_type": "E",
            "thermo_kwargs": {"temperature": 298.15},
            "threshold": None,
            "matrix_format": "xlsx",
        }
        create_grouper_job_from_context(ctx, strategy="rmsd")
        _, call_kwargs = mock_job_cls.call_args
        assert call_kwargs["thermo_parameters"] is None
