from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob


def _invoke_xtb_and_capture_settings(job_path, args):
    runner = CliRunner()
    mock_job = MagicMock()
    with patch(job_path, return_value=mock_job) as mock_job_cls:
        result = runner.invoke(
            xtb,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=False,
        )
        settings = mock_job_cls.call_args[1]["settings"]
    return result, settings, mock_job_cls


class TestXTBHelp:
    def test_run_help_lists_xtb(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "xtb" in result.output

    def test_sub_help_lists_xtb(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "xtb" in result.output

    def test_xtb_help_lists_jobs(self):
        result = CliRunner().invoke(run, ["xtb", "--help"])
        assert result.exit_code == 0, result.output
        assert "opt" in result.output
        assert "sp" in result.output
        assert "hess" in result.output


class TestXTBCLISettings:
    def test_sp_settings_merge(self, single_molecule_xyz_file):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "-1",
                "-m",
                "2",
                "-g",
                "gfn1",
                "-sm",
                "alpb",
                "-si",
                "water",
                "--grad",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_args[0] == ()
        assert settings.jobtype == "sp"
        assert settings.charge == -1
        assert settings.multiplicity == 2
        assert settings.gfn_version == "gfn1"
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "water"
        assert settings.grad is True
        assert settings.freq is False

    def test_opt_instantiates_xtb_opt_job(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.opt.XTBOptJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
                "-O",
                "loose",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "opt"
        assert settings.optimization_level == "loose"
        assert settings.freq is True

    def test_opt_constrain_and_force_constant(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.opt.XTBOptJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
                "-c",
                "[[1,2],[1,2,3],[1,2,3,4]]",
                "-f",
                "0.25",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "opt"
        assert settings.constraints == [[1, 2], [1, 2, 3], [1, 2, 3, 4]]
        assert settings.force_constant == 0.25

    def test_hess_instantiates_xtb_hess_job(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.hess.XTBHessJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "hess",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "hess"

    def test_sp_additional_flags_merge(self, single_molecule_xyz_file):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "--copy --json",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_flags == "--copy --json"

    def test_settings_from_xtb_main_output(self, xtb_p_benzyne_opt_outfolder):
        filepath = (
            f"{xtb_p_benzyne_opt_outfolder}/" f"p_benzyne_opt_alpb_toluene.out"
        )
        settings = XTBJobSettings.from_filepath(filepath)

        assert settings.gfn_version == "gfn2"
        assert settings.optimization_level == "loose"
        assert settings.charge == 0
        assert settings.multiplicity == 3
        assert settings.jobtype == "opt"
        assert settings.freq is False
        assert settings.grad is True
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "toluene"

    def test_settings_from_xtb_sp_output(self, xtb_p_benzyne_sp_outfolder):
        filepath = (
            f"{xtb_p_benzyne_sp_outfolder}/" f"p_benzyne_sp_alpb_toluene.out"
        )
        settings = XTBJobSettings.from_filepath(filepath)

        assert settings.gfn_version == "gfn2"
        assert settings.optimization_level is None
        assert settings.charge == 0
        assert settings.multiplicity == 3
        assert settings.jobtype == "sp"
        assert settings.freq is False
        assert settings.grad is False
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "toluene"

    def test_hess_accepts_xtb_main_output(self, xtb_co2_outfolder):
        xtb_main_output = str(Path(xtb_co2_outfolder) / "co2_ohess.out")
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.hess.XTBHessJob",
            ["-p", "test", "-f", xtb_main_output, "hess"],
        )

        assert result.exit_code == 0, result.output
        assert settings.jobtype == "hess"
        assert settings.freq is True
        assert settings.charge == 0
        assert settings.multiplicity == 1
        assert (
            mock_job_cls.call_args.kwargs["molecule"].empirical_formula
            == "CO2"
        )
        assert mock_job_cls.call_args.kwargs["label"] == "co2_ohess_hess"


class TestXTBCLIGroupValidation:
    """Validation and settings-merge branches on the ``xtb`` group
    callback itself (not delegated to any subcommand)."""

    def test_molecule_id_not_supported_raises(self, single_molecule_xyz_file):
        result = CliRunner().invoke(
            xtb,
            [
                "-f",
                single_molecule_xyz_file,
                "--mid",
                "abc",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "not supported for xTB job submission" in result.output

    def test_index_and_structure_index_mutually_exclusive(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            xtb,
            [
                "-f",
                single_molecule_xyz_file,
                "-i",
                "1",
                "--si",
                "1",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_structure_index_used_as_index(self, multiple_molecules_xyz_file):
        """``--si`` is treated as an alias for ``-i`` when given alone."""
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                multiple_molecules_xyz_file,
                "--si",
                "1",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_requires_exactly_one_selector(
        self, database_chemsmart_file
    ):
        result = CliRunner().invoke(
            xtb,
            ["-f", database_chemsmart_file, "sp"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_rejects_multiple_selectors(
        self, database_chemsmart_file
    ):
        result = CliRunner().invoke(
            xtb,
            [
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "--sid",
                "abc",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_index_requires_record_selector(
        self, database_chemsmart_file
    ):
        result = CliRunner().invoke(
            xtb,
            [
                "-f",
                database_chemsmart_file,
                "--sid",
                "abc",
                "-i",
                "1",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "can only be used together with" in result.output

    def test_non_chemsmart_db_falls_back_to_defaults(self, database_ase_file):
        """A .db file that isn't a chemsmart database (e.g. a plain ASE
        db) falls back to default xTB settings entirely."""
        mock_molecule = MagicMock(name="ase_db_molecule")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=[mock_molecule],
        ):
            result, settings, _ = _invoke_xtb_and_capture_settings(
                "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
                [
                    "-p",
                    "test",
                    "-f",
                    database_ase_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "sp",
                ],
            )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_filename_and_pubchem_both_missing_raises(self):
        result = CliRunner().invoke(
            xtb,
            ["-p", "test", "sp"],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "has not been specified" in str(result.exception)

    def test_filename_and_pubchem_both_given_raises(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            xtb,
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "--pubchem",
                "222",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "have been specified" in str(result.exception)

    def test_label_and_append_label_mutually_exclusive_raises(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            xtb,
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-l",
                "custom",
                "-a",
                "suffix",
                "sp",
            ],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "not both" in str(result.exception)

    def test_append_label_suffixes_filename_derived_label(
        self, single_molecule_xyz_file
    ):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-a",
                "suffix",
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label.endswith("_suffix")

    def test_pubchem_only_without_label_falls_back_to_output(self):
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with (
            patch(
                "chemsmart.io.molecules.structure.Molecule.from_pubchem",
                return_value=[pubchem_molecule],
            ),
            patch(
                "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob"
            ) as mock_job_cls,
        ):
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                xtb,
                [
                    "-p",
                    "test",
                    "--pubchem",
                    "222",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "sp",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label == "output_sp"

    def test_pubchem_only_with_append_label_falls_back_to_output(self):
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with (
            patch(
                "chemsmart.io.molecules.structure.Molecule.from_pubchem",
                return_value=[pubchem_molecule],
            ),
            patch(
                "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob"
            ) as mock_job_cls,
        ):
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                xtb,
                [
                    "-p",
                    "test",
                    "--pubchem",
                    "222",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "-a",
                    "suffix",
                    "sp",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label == "output_suffix"

    def test_default_label_suffixed_with_invoked_subcommand(
        self, single_molecule_xyz_file
    ):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
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
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label.endswith("_sp")

    def test_chemsmart_db_loads_molecule_by_record_index(
        self, database_chemsmart_file
    ):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_loads_molecule_by_structure_id(
        self, database_chemsmart_file
    ):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_append_label_includes_structure_id_suffix(
        self, database_chemsmart_file
    ):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "-a",
                "suffix",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "SID-f751bb2c27e2" in label

    def test_chemsmart_db_append_label_includes_record_id_suffix(
        self, database_chemsmart_file
    ):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--rid",
                "6de213a0",
                "-a",
                "suffix",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "RID-6de213a0" in label

    def test_chemsmart_db_default_label_includes_record_id_suffix(
        self, database_chemsmart_file
    ):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--rid",
                "6de213a0",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "RID-6de213a0" in label

    def test_chemsmart_db_append_label_includes_record_index_suffix(
        self, database_chemsmart_file
    ):
        result, settings, mock_job_cls = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "-a",
                "suffix",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "RI-1" in label

    def test_remove_solvent_clears_solvent_fields(
        self, single_molecule_xyz_file
    ):
        result, settings, _ = _invoke_xtb_and_capture_settings(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "gbsa",
                "-si",
                "water",
                "--remove-solvent",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self, multiple_molecules_xyz_file
    ):
        with patch(
            "chemsmart.jobs.xtb.singlepoint.XTBSinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                xtb,
                [
                    "-p",
                    "test",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "sp",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]


class TestXTBJobClasses:
    def test_xtb_job_types_are_distinct(self):
        assert XTBOptJob.TYPE == "xtbopt"
        assert XTBSinglePointJob.TYPE == "xtbsp"
        assert XTBHessJob.TYPE == "xtbhess"


class TestXTBSubmission:
    def test_sub_test_writes_xtb_scripts(
        self, tmp_path, server_yaml_file, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        xyz_file = (
            Path(__file__).resolve().parent
            / "data"
            / "XTBTests"
            / "inputs"
            / "water.xyz"
        )
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "xtb",
                "-p",
                "test",
                "-f",
                str(xyz_file),
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output
        run_script = tmp_path / "chemsmart_run_water_sp.py"
        submit_script = tmp_path / "chemsmart_sub_water_sp.sh"
        assert run_script.exists()
        assert submit_script.exists()
        assert "'--test'" not in run_script.read_text()
        assert "'xtb'" in run_script.read_text()
        assert (
            "conda activate ~/anaconda3/envs/chemsmart"
            in submit_script.read_text()
        )
