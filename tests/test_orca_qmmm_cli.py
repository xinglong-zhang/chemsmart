"""
Direct tests for the ``qmmm`` subcommand created by
``create_orca_qmmm_subcommand`` in ``chemsmart.cli.orca.qmmm``.

The subcommand is attached to several ORCA jobtype groups (``opt``,
``ts``, ``sp``, ``scan``, ``qrc``, ``modred``, ``neb``); ``opt`` is
used as the parent throughout since it's the simplest.
"""

from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from chemsmart.cli.orca.orca import orca


class TestOrcaQmmmSubcommand:
    def test_all_qmmm_options_applied_to_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
                "qmmm",
                "-j",
                "QM/QM2/MM",
                "-hx",
                "b3lyp",
                "-hb",
                "def2-svp",
                "-ix",
                "pbe",
                "-ib",
                "def2-svp",
                "-im",
                "XTB",
                "-lm",
                "AMBER",
                "-ha",
                "1-3",
                "-ia",
                "4-5",
                "-ct",
                "0",
                "-mt",
                "1",
                "-ci",
                "0",
                "-mi",
                "1",
                "-ch",
                "0",
                "-mh",
                "1",
                "-s",
                "CPCM",
                "-a",
                "1-3",
                "-ua",
                "pdb_info",
                "-o",
                "1-2",
                "-d",
                "true",
                "-db",
                "true",
                "-e",
                "electronic",
                "-cc",
                "true",
                "-xn",
                "50",
                "-t",
                "0.0001",
                "-sc",
                "0.5",
                "-nc",
                "10",
                "-ecp",
                "def2-ecp",
                "-ecpn",
                "2",
                "-sc2",
                "0.5",
            ],
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.jobtype == "QM/QM2/MM"
        assert settings.high_level_functional == "b3lyp"
        assert settings.high_level_basis == "def2-svp"
        assert settings.intermediate_level_functional == "pbe"
        assert settings.intermediate_level_basis == "def2-svp"
        assert settings.intermediate_level_method == "XTB"
        assert settings.low_level_method == "AMBER"
        assert settings.charge_total == 0
        assert settings.mult_total == 1
        assert settings.charge_intermediate == 0
        assert settings.mult_intermediate == 1
        assert settings.charge_high == 0
        assert settings.mult_high == 1
        assert settings.intermediate_level_solvation == "CPCM"
        assert settings.optregion_fixed_atoms == "1-2"
        assert settings.use_active_info_from_pbc == "pdb_info"
        assert settings.delete_la_double_counting is True
        assert settings.delete_la_bond_double_counting_atoms is True
        assert settings.embedding_type == "electronic"
        assert settings.conv_charges is True
        assert settings.conv_charges_max_n_cycles == 50
        assert settings.conv_charges_conv_thresh == 0.0001
        assert settings.scale_formal_charge_mm_atom == 0.5
        assert settings.n_unit_cell_atoms == 10
        assert settings.ecp_layer_ecp == "def2-ecp"
        assert settings.ecp_layer == 2
        assert settings.scale_formal_charge_ecp_atom == 0.5

    def test_project_qmmm_settings_loaded_from_yaml(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """The 'test_qmmm' project's qmmm: section should seed the
        settings (non-None qmmm_settings() branch)."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "test_qmmm",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
                "qmmm",
                "-ha",
                "1-3",
            ],
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_charge_and_multiplicity_populated_from_intermediate(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """_populate_charge_and_multiplicity_on_settings should copy
        charge_intermediate/mult_intermediate onto charge/multiplicity
        when they're the most specific pair given."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-ci",
                "-1",
                "-mi",
                "2",
            ],
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings.charge == -1
        assert settings.multiplicity == 2

    def test_charge_and_multiplicity_populated_from_high(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-ch",
                "1",
                "-mh",
                "3",
            ],
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings.charge == 1
        assert settings.multiplicity == 3

    def test_charge_and_multiplicity_populated_from_total(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-ct",
                "2",
                "-mt",
                "1",
            ],
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings.charge == 2
        assert settings.multiplicity == 1

    def test_high_level_and_intermediate_atoms_attached_to_molecule(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        from chemsmart.jobs.orca.qmmm import ORCAQMMMJob

        with patch.object(
            ORCAQMMMJob, "__init__", return_value=None
        ) as mock_init:
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "opt",
                    "qmmm",
                    "-ha",
                    "1-2",
                    "-ia",
                    "3",
                ],
                obj={"jobrunner": None},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_init.call_count == 1
        _, kwargs = mock_init.call_args
        molecule = kwargs["molecule"]
        assert molecule.high_level_atoms is not None
        assert molecule.intermediate_level_atoms is not None

    def test_label_already_containing_qmmm_not_double_suffixed(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        from chemsmart.jobs.orca.qmmm import ORCAQMMMJob

        with patch.object(
            ORCAQMMMJob, "__new__", return_value=MagicMock()
        ) as mock_new:
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-l",
                    "myjob_qmmm",
                    "opt",
                    "qmmm",
                    "-ha",
                    "1-3",
                ],
                obj={"jobrunner": None},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        _, kwargs = mock_new.call_args
        assert kwargs["label"] == "myjob_qmmm"

    def test_high_level_h_bond_length_option_is_unusable(
        self,
        single_molecule_xyz_file,
    ):
        """BUG (see BUGS_FOUND.md #28): -h/--high-level-h-bond-length is
        declared with click's type=dict, which converts the raw CLI
        string via Python's dict(value) constructor -- not a literal
        dict parser. Any non-empty value therefore fails at the click
        parsing layer itself, before ast.literal_eval (used later on
        this same value) is ever reached, making the option unusable
        for its documented purpose."""
        result = CliRunner().invoke(
            orca,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-h",
                "{1: 1.1}",
            ],
            obj={"jobrunner": None},
        )
        assert result.exit_code != 0
        assert "Invalid value for '-h'" in result.output
