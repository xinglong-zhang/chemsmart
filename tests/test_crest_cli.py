import os
from pathlib import Path
from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.crest.crest import crest
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.crest import CRESTConformerSearchJob


class TestCRESTHelp:
    def test_run_help_lists_crest(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "crest" in result.output

    def test_sub_help_lists_crest(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "crest" in result.output

    def test_crest_help_lists_subcommands(self):
        result = CliRunner().invoke(run, ["crest", "--help"])
        assert result.exit_code == 0, result.output
        assert "conformers" in result.output


class TestCRESTCLISettings:
    def test_conformers_settings_merge(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
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
                "-O",
                "loose",
                "-sm",
                "alpb",
                "-si",
                "water",
                "-w",
                "8.0",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "conformers"
        assert settings.charge == -1
        assert settings.multiplicity == 2
        assert settings.gfn_version == "gfn1"
        assert settings.optimization_level == "loose"
        assert settings.solvent_model == "alpb"
        assert settings.solvent_id == "water"
        assert settings.energy_window == 8.0

    def test_conformers_default_settings(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.jobtype == "conformers"
        assert settings.gfn_version == "gfn2"
        assert settings.constraints is None

    def test_constrained_settings(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
                "-c",
                "[[1,2],[1,3]]",
                "--force-constant",
                "0.5",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.constraints == [[1, 2], [1, 3]]
        assert settings.force_constant == 0.5


class TestCRESTCLILabelFallback:
    """Covers the filename-less label fallback branches in the ``crest``
    group callback (no -f/--filename, molecule obtained via --pubchem)."""

    def test_pubchem_only_without_label_falls_back_to_crest_output(self):
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with (
            patch(
                "chemsmart.io.molecules.structure.Molecule.from_pubchem",
                return_value=[pubchem_molecule],
            ),
            patch(
                "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob"
            ) as mock_job_cls,
        ):
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                crest,
                [
                    "-p",
                    "test",
                    "--pubchem",
                    "222",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "conformers",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label == "crest_output_conformers"

    def test_pubchem_only_with_append_label_falls_back_to_output(self):
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with (
            patch(
                "chemsmart.io.molecules.structure.Molecule.from_pubchem",
                return_value=[pubchem_molecule],
            ),
            patch(
                "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob"
            ) as mock_job_cls,
        ):
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                crest,
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
                    "conformers",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label == "output_suffix"

    def test_append_label_suffixes_filename_derived_label(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, _, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
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
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label.endswith("_suffix")


class TestCRESTJobClasses:
    def test_crest_job_types_are_distinct(self):
        assert CRESTConformerSearchJob.TYPE == "crestconformers"


class TestCRESTCLIGroupValidation:
    """Validation and settings-merge branches on the ``crest`` group
    callback itself (not delegated to the ``conformers`` subcommand)."""

    def test_molecule_id_not_supported_raises(self, single_molecule_xyz_file):
        result = CliRunner().invoke(
            crest,
            [
                "-f",
                single_molecule_xyz_file,
                "--mid",
                "abc",
                "conformers",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "not supported for CREST job submission" in result.output

    def test_index_and_structure_index_mutually_exclusive(
        self, single_molecule_xyz_file
    ):
        result = CliRunner().invoke(
            crest,
            [
                "-f",
                single_molecule_xyz_file,
                "-i",
                "1",
                "--si",
                "1",
                "conformers",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_structure_index_used_as_index(
        self, multiple_molecules_xyz_file, run_crest_and_capture_settings
    ):
        """``--si`` is treated as an alias for ``-i`` when given alone."""
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
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
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_requires_exactly_one_selector(
        self, database_chemsmart_file
    ):
        result = CliRunner().invoke(
            crest,
            ["-f", database_chemsmart_file, "conformers"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_rejects_multiple_selectors(
        self, database_chemsmart_file
    ):
        result = CliRunner().invoke(
            crest,
            [
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "--sid",
                "abc",
                "conformers",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_index_requires_record_selector(
        self, database_chemsmart_file
    ):
        result = CliRunner().invoke(
            crest,
            [
                "-f",
                database_chemsmart_file,
                "--sid",
                "abc",
                "-i",
                "1",
                "conformers",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code != 0
        assert "can only be used together with" in result.output

    def test_xtb_output_inherits_charge_and_multiplicity(
        self, xtb_water_outfolder, run_crest_and_capture_settings
    ):
        """A ``.out`` file detected as xTB should fall back to default
        CREST settings but inherit charge/multiplicity from xTB."""
        xtb_out = os.path.join(xtb_water_outfolder, "water_ohess.out")
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            ["-p", "test", "-f", xtb_out, "conformers"],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_non_chemsmart_db_falls_back_to_defaults(
        self, database_ase_file, run_crest_and_capture_settings
    ):
        """A .db file that isn't a chemsmart database (e.g. a plain ASE
        db) falls back to default CREST settings entirely."""
        mock_molecule = MagicMock(name="ase_db_molecule")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=[mock_molecule],
        ):
            result, settings, _ = run_crest_and_capture_settings(
                "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
                [
                    "-p",
                    "test",
                    "-f",
                    database_ase_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "conformers",
                ],
            )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_filename_and_pubchem_both_missing_raises(self):
        result = CliRunner().invoke(
            crest,
            ["-p", "test", "conformers"],
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
            crest,
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "--pubchem",
                "222",
                "conformers",
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
            crest,
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-l",
                "custom",
                "-a",
                "suffix",
                "conformers",
            ],
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "not both" in str(result.exception)

    def test_append_label_suffixes_filename_derived_label(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
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
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label.endswith("_suffix")

    def test_pubchem_only_with_append_label_falls_back_to_output(self):
        """Previously a bug native to crest.py (not merge-related): the
        append_label branch called os.path.basename(filename)
        unconditionally, without the `if filename:` guard that the
        default-label branch just below it already has, crashing with a
        TypeError for a filename-less (PubChem-only) job combined with
        -a/--append-label. Fixed by adding the same guard."""
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with (
            patch(
                "chemsmart.io.molecules.structure.Molecule.from_pubchem",
                return_value=[pubchem_molecule],
            ),
            patch(
                "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob"
            ) as mock_job_cls,
        ):
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                crest,
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
                    "conformers",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label == "output_suffix"

    def test_default_label_suffixed_with_invoked_subcommand(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label.endswith("_conformers")

    def test_chemsmart_db_loads_molecule_by_record_index(
        self, database_chemsmart_file, run_crest_and_capture_settings
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_loads_molecule_by_structure_id(
        self, database_chemsmart_file, run_crest_and_capture_settings
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_append_label_includes_structure_id_suffix(
        self, database_chemsmart_file, run_crest_and_capture_settings
    ):
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "-a",
                "suffix",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "SID-f751bb2c27e2" in label

    def test_chemsmart_db_append_label_includes_record_index_suffix(
        self, database_chemsmart_file, run_crest_and_capture_settings
    ):
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "-a",
                "suffix",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "RI-1" in label

    def test_chemsmart_db_append_label_includes_record_id_suffix(
        self, database_chemsmart_file, run_crest_and_capture_settings
    ):
        """--rid (record_id, no structure_id) takes the elif branch for
        the append-label suffix."""
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--rid",
                "6de213a0",
                "-a",
                "suffix",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "RID-6de213a0" in label

    def test_chemsmart_db_default_label_includes_record_id_suffix(
        self, database_chemsmart_file, run_crest_and_capture_settings
    ):
        """--rid alone (no -a/-l) takes the default-label elif branch
        for the record_id suffix."""
        result, settings, mock_job_cls = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                database_chemsmart_file,
                "--rid",
                "6de213a0",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert "RID-6de213a0" in label

    def test_nci_flag_enabled(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--nci",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.nci is True

    def test_additional_flags_option(
        self, single_molecule_xyz_file, run_crest_and_capture_settings
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
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
                "--squick",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_flags == "--squick"

    def test_remove_solvent_clears_solvent_fields(
        self,
        crest_yaml_settings_solv_project_name,
        single_molecule_xyz_file,
        run_crest_and_capture_settings,
    ):
        result, settings, _ = run_crest_and_capture_settings(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            [
                "-p",
                crest_yaml_settings_solv_project_name,
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "conformers",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self, multiple_molecules_xyz_file
    ):
        with patch(
            "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                crest,
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
                    "conformers",
                ],
                obj={"jobrunner": MagicMock()},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]


class TestCRESTSubmission:
    def test_sub_test_writes_crest_scripts(
        self, tmp_path, server_yaml_file, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        xyz_file = (
            Path(__file__).resolve().parent
            / "data"
            / "StructuresTests"
            / "xyz"
            / "1-mer.xyz"
        )
        result = CliRunner().invoke(
            sub,
            [
                "-s",
                server_yaml_file,
                "--test",
                "crest",
                "-p",
                "test",
                "-f",
                str(xyz_file),
                "-c",
                "0",
                "-m",
                "1",
                "conformers",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, result.output
        run_script = tmp_path / "chemsmart_run_crest_1-mer_conformers.py"
        submit_script = tmp_path / "chemsmart_sub_crest_1-mer_conformers.sh"
        assert run_script.exists()
        assert submit_script.exists()
        assert "'--test'" not in run_script.read_text()
        assert "'crest'" in run_script.read_text()
        assert (
            "conda activate ~/anaconda3/envs/chemsmart"
            in submit_script.read_text()
        )
