import importlib
import logging
import os
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.thermochemistry.thermochemistry import thermochemistry
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.runner import FakeGaussianJobRunner
from chemsmart.jobs.nciplot.runner import FakeNCIPLOTJobRunner
from chemsmart.jobs.orca.runner import FakeORCAJobRunner
from chemsmart.settings.server import Server

thermochemistry_cli_module = importlib.import_module(
    "chemsmart.cli.thermochemistry.thermochemistry"
)

mol_cli_module = importlib.import_module("chemsmart.cli.mol.mol")


############ IO Fixtures ####################################


############ Thermochemistry Mock Fixtures ##################
@pytest.fixture()
def make_thermochemistry_mock():
    """Factory fixture that creates a MagicMock mimicking a Thermochemistry instance.

    Returns a callable that accepts the attributes accessed by
    ``Thermochemistry.cleaned_frequencies`` and returns a properly configured
    mock, so test methods stay free of direct ``MagicMock`` construction.

    Usage::

        def test_something(make_thermochemistry_mock):
            mock = make_thermochemistry_mock(
                vibrational_frequencies=[-50.0, 100.0],
                jobtype="opt",
                check_imaginary_frequencies=False,
            )
            result = Thermochemistry.cleaned_frequencies.fget(mock)
            assert result == [100.0, 100.0]
    """
    from chemsmart.analysis.thermochemistry import (
        NEAR_ZERO_FREQUENCY_TOLERANCE_CM,
        Thermochemistry,
    )

    def _factory(
        vibrational_frequencies,
        jobtype="opt",
        check_imaginary_frequencies=True,
        s_freq_cutoff_cm=None,
        h_freq_cutoff_cm=None,
        rotational_mode="gaussian",
        near_zero_frequency_tolerance_cm=NEAR_ZERO_FREQUENCY_TOLERANCE_CM,
        reaction_coordinate_mode=0,
    ):
        mock = MagicMock(spec=Thermochemistry)
        mock.vibrational_frequencies = vibrational_frequencies
        mock.jobtype = jobtype
        mock.check_imaginary_frequencies = check_imaginary_frequencies
        mock.s_freq_cutoff_cm = s_freq_cutoff_cm
        mock.h_freq_cutoff_cm = h_freq_cutoff_cm
        mock.rotational_mode = rotational_mode
        # Instance attributes set in __init__ are invisible to spec=; the
        # near-zero convention's tolerance must be mirrored explicitly.
        mock.near_zero_frequency_tolerance_cm = (
            near_zero_frequency_tolerance_cm
        )
        # Which printed mode the session named as the reaction
        # coordinate, 1-based; 0 means the first genuine imaginary one,
        # which is what every permissive expectation below pins.
        mock.reaction_coordinate_mode = reaction_coordinate_mode
        mock.filename = "dummy.log"
        mock.target = "dummy.log"
        return mock

    return _factory


############ CLI Fixtures ##################
@pytest.fixture
def make_cli_ctx_obj():
    """Factory for the minimal Click context object."""

    def _make(jobrunner):
        return {"jobrunner": jobrunner}

    return _make


@pytest.fixture()
def invoke_config_server():
    """Return a callable that invokes 'chemsmart config server' via Click's CliRunner.

    Usage in tests::

        def test_something(invoke_config_server):
            result = invoke_config_server()
            assert result.exit_code == 0
    """
    from chemsmart.cli.config import config

    def _invoke(args=None):
        runner = CliRunner()
        return runner.invoke(config, ["server"] + (args or []))

    return _invoke


@pytest.fixture()
def run_thermochemistry_and_capture_settings():
    """Run the thermochemistry CLI with mocked job construction."""

    def _run(
        extra_args=None,
        ctx_obj=None,
        filename="dummy.log",
        detected_program="gaussian",
    ):
        runner = CliRunner()
        captured_settings = None
        mock_job = MagicMock()

        base_args = ["-f", filename, "-T", "298.15"]
        cli_args = base_args + (extra_args or [])

        with (
            patch.object(
                thermochemistry_cli_module,
                "get_program_type_from_file",
                return_value=detected_program,
            ),
            patch.object(
                thermochemistry_cli_module.ThermochemistryJob,
                "from_filename",
                return_value=mock_job,
            ) as mock_from_filename,
        ):
            result = runner.invoke(
                thermochemistry,
                cli_args,
                obj=ctx_obj or {},
                catch_exceptions=False,
            )
            if mock_from_filename.call_args is not None:
                captured_settings = mock_from_filename.call_args[1].get(
                    "settings"
                )

        return result, captured_settings

    return _run


@pytest.fixture()
def run_thermochemistry_with_directory():
    """Fixture to invoke thermochemistry CLI with directory options and mocked folder.

    Patches ``BaseFolder`` so that
    ``get_all_output_files_in_current_folder_by_program``,
    ``get_all_files_in_current_folder_by_suffix``, and
    ``get_all_files_in_current_folder_by_program_and_suffix`` all return
    the caller-supplied ``mock_files`` list.  Also patches
    ``ThermochemistryJob.from_filename`` to avoid real job execution.

    Usage::

        def test_something(run_thermochemistry_with_directory, tmp_path):
            result, mock_from_filename = run_thermochemistry_with_directory(
                ["-d", str(tmp_path), "-p", "gaussian", "-T", "298.15"],
                mock_files=["/fake/a.log", "/fake/b.log"],
            )
            assert result.exit_code == 0
            assert mock_from_filename.call_count == 2
    """

    def _invoke(extra_args, mock_files=None):
        if mock_files is None:
            mock_files = ["/fake/dir/mol1.log"]

        runner = CliRunner()
        mock_job = MagicMock()
        mock_job.label = "mol1"

        with (
            patch.object(
                thermochemistry_cli_module, "BaseFolder"
            ) as mock_folder_cls,
            patch.object(
                thermochemistry_cli_module.ThermochemistryJob, "from_filename"
            ) as mock_from_filename,
        ):
            mock_folder = MagicMock()
            # Configure every discovery method to return the caller-supplied list
            mock_folder.get_all_output_files_in_current_folder_by_program.return_value = (
                mock_files
            )
            mock_folder.get_all_output_files_in_current_folder_and_subfolders_by_program.return_value = (
                mock_files
            )
            mock_folder.get_all_files_in_current_folder_by_suffix.return_value = (
                mock_files
            )
            mock_folder.get_all_files_in_current_folder_by_program_and_suffix.return_value = (
                mock_files
            )
            mock_folder_cls.return_value = mock_folder
            mock_from_filename.return_value = mock_job

            result = runner.invoke(thermochemistry, extra_args)
            return result, mock_from_filename

    return _invoke


@pytest.fixture()
def run_gaussian_and_capture_settings():
    """Run the gaussian CLI with a patched job class and capture settings."""

    def _run(job_class_path, cli_args, ctx_obj):
        runner = CliRunner()
        captured_settings = None

        with patch(job_class_path) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
                cli_args,
                obj=ctx_obj,
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured_settings = mock_job_cls.call_args[1].get("settings")

        return result, captured_settings

    return _run


@pytest.fixture()
def run_orca_and_capture_settings():
    """Run the orca CLI with a patched job class and capture settings."""
    from chemsmart.cli.orca.orca import orca as orca_cli

    def _run(job_class_path, cli_args, ctx_obj=None):
        if ctx_obj is None:
            ctx_obj = {}
        runner = CliRunner()
        captured_settings = None

        with patch(job_class_path) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                cli_args,
                obj=ctx_obj,
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured_settings = mock_job_cls.call_args[1].get("settings")

        return result, captured_settings

    return _run


############ Gaussian Fixtures ##################
@pytest.fixture()
def test_data_directory():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(current_directory, "data"))


# master gaussian test directory
@pytest.fixture()
def gaussian_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "GaussianTests")


# Gaussian folder for semiempirical calculations
@pytest.fixture()
def gaussian_semiempirical_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "semiempirical")


@pytest.fixture()
def gaussian_semiempirical_pm6_output_file(
    gaussian_semiempirical_test_directory,
):
    return os.path.join(
        gaussian_semiempirical_test_directory, "DBU_PM6_opt.log"
    )


# Gaussian output file from outputs folder
@pytest.fixture()
def outputs_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "outputs")


@pytest.fixture()
def wbi_outputfile(outputs_test_directory):
    wbi_outputfile = os.path.join(
        outputs_test_directory, "TS_5coord_XIII_wbi.log"
    )
    return wbi_outputfile


# Gaussian input files
@pytest.fixture()
def gaussian_inputs_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "inputs")


@pytest.fixture()
def gaussian_opt_inputfile(gaussian_inputs_test_directory):
    gaussian_opt_input = os.path.join(
        gaussian_inputs_test_directory, "model_opt_input.com"
    )
    return gaussian_opt_input


@pytest.fixture()
def gaussian_frozen_opt_inputfile(gaussian_inputs_test_directory):
    gaussian_frozen_opt_inputfile = os.path.join(
        gaussian_inputs_test_directory, "frozen_coordinates_opt.com"
    )
    return gaussian_frozen_opt_inputfile


@pytest.fixture()
def gaussian_modred_inputfile(gaussian_inputs_test_directory):
    gaussian_modred_inputfile = os.path.join(
        gaussian_inputs_test_directory, "model_modred_input.com"
    )
    return gaussian_modred_inputfile


@pytest.fixture()
def gaussian_scan_inputfile(gaussian_inputs_test_directory):
    gaussian_scan_inputfile = os.path.join(
        gaussian_inputs_test_directory, "model_scan_input.com"
    )
    return gaussian_scan_inputfile


@pytest.fixture()
def hf_com_filepath(gaussian_inputs_test_directory):
    return os.path.join(gaussian_inputs_test_directory, "hf.com")


# Gaussian input files for genecp
@pytest.fixture()
def gaussian_inputs_genecp_directory(gaussian_inputs_test_directory):
    return os.path.join(gaussian_inputs_test_directory, "genecp")


@pytest.fixture()
def gaussian_opt_genecp_inputfile(gaussian_inputs_genecp_directory):
    gaussian_opt_genecp_input = os.path.join(
        gaussian_inputs_genecp_directory, "opt_genecp.com"
    )
    return gaussian_opt_genecp_input


@pytest.fixture()
def modred_gen_inputfile(gaussian_inputs_genecp_directory):
    return os.path.join(gaussian_inputs_genecp_directory, "modred_gen.com")


@pytest.fixture()
def modred_genecp_inputfile(gaussian_inputs_genecp_directory):
    return os.path.join(gaussian_inputs_genecp_directory, "modred_genecp.com")


@pytest.fixture()
def modred_genecp_custom_solvent_inputfile(gaussian_inputs_genecp_directory):
    return os.path.join(
        gaussian_inputs_genecp_directory, "modred_genecp_custom_solvent.com"
    )


# Gaussian input files for link jobs
@pytest.fixture()
def gaussian_link_inputs_test_directory(gaussian_inputs_test_directory):
    return os.path.join(gaussian_inputs_test_directory, "link")


@pytest.fixture()
def gaussian_link_opt_input(gaussian_link_inputs_test_directory):
    return os.path.join(
        gaussian_link_inputs_test_directory, "link_opt_input_opt_link.com"
    )


@pytest.fixture()
def gaussian_link_ts_input(gaussian_link_inputs_test_directory):
    return os.path.join(
        gaussian_link_inputs_test_directory, "link_ts_input_ts_link.com"
    )


# Gaussian output files for link jobs
@pytest.fixture()
def gaussian_link_outputs_test_directory(gaussian_outputs_test_directory):
    gaussian_link_outputs_test_directory = os.path.join(
        gaussian_outputs_test_directory, "link"
    )
    return gaussian_link_outputs_test_directory


@pytest.fixture()
def gaussian_link_opt_outputfile(gaussian_link_outputs_test_directory):
    gaussian_link_opt_outfile = os.path.join(
        gaussian_link_outputs_test_directory,
        "oxygen_openshell_singlet_opt_link.log",
    )
    return gaussian_link_opt_outfile


@pytest.fixture()
def gaussian_link_ts_outputfile(gaussian_link_outputs_test_directory):
    gaussian_link_ts_outfile = os.path.join(
        gaussian_link_outputs_test_directory,
        "oxygen_openshell_singlet_ts_link.log",
    )
    return gaussian_link_ts_outfile


@pytest.fixture()
def gaussian_link_modred_output(gaussian_link_outputs_test_directory):
    gaussian_link_modred_outfile = os.path.join(
        gaussian_link_outputs_test_directory,
        "fe_ch_quintet_modred_link.log",
    )
    return gaussian_link_modred_outfile


@pytest.fixture()
def gaussian_link_sp_outputfile(gaussian_link_outputs_test_directory):
    return os.path.join(
        gaussian_link_outputs_test_directory,
        "oxygen_openshell_singlet_sp_link.log",
    )


@pytest.fixture()
def gaussian_failed_link_output(gaussian_link_outputs_test_directory):
    return os.path.join(
        gaussian_link_outputs_test_directory, "failed_link_job.log"
    )


@pytest.fixture()
def gaussian_link_sp_input(gaussian_link_inputs_test_directory):
    return os.path.join(
        gaussian_link_inputs_test_directory, "link_sp_input_sp_link.com"
    )


@pytest.fixture()
def gaussian_qmmm_input_test_directory(gaussian_inputs_test_directory):
    return os.path.join(gaussian_inputs_test_directory, "qmmm")


@pytest.fixture()
def gaussian_qmmm_inputfile_2layer(gaussian_qmmm_input_test_directory):
    return os.path.join(gaussian_qmmm_input_test_directory, "CH3CH3.com")


@pytest.fixture()
def gaussian_qmmm_inputfile_3layer(gaussian_qmmm_input_test_directory):
    return os.path.join(gaussian_qmmm_input_test_directory, "CH3COOH.com")


# Gaussian output files
@pytest.fixture()
def gaussian_outputs_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "outputs")


@pytest.fixture()
def gaussian_singlet_opt_outfile(gaussian_outputs_test_directory):
    gaussian_singlet_opt_output = os.path.join(
        gaussian_outputs_test_directory, "nhc_neutral_singlet.log"
    )
    return gaussian_singlet_opt_output


@pytest.fixture()
def gaussian_triplet_opt_outfile(gaussian_outputs_test_directory):
    gaussian_triplet_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "iron_neutral_triplet.log"
    )
    return gaussian_triplet_opt_outfile


@pytest.fixture()
def gaussian_quintet_opt_outfile(gaussian_outputs_test_directory):
    gaussian_quintet_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "iron_neutral_quintet.log"
    )
    return gaussian_quintet_opt_outfile


# Gaussian output files for genecp
@pytest.fixture()
def gaussian_ts_genecp_outfile(gaussian_outputs_test_directory):
    gaussian_ts_genecp_output = os.path.join(
        gaussian_outputs_test_directory, "pd_genecp_ts.log"
    )
    return gaussian_ts_genecp_output


@pytest.fixture()
def gaussian_pd_insertion_ts_r_outfile(gaussian_outputs_test_directory):
    return os.path.join(
        gaussian_outputs_test_directory, "Pd_insertion_ts_r.log"
    )


@pytest.fixture()
def gaussian_full_gen_outfile(gaussian_outputs_test_directory):
    return os.path.join(
        gaussian_outputs_test_directory,
        "bromochloromethane_full_gen.log",
    )


@pytest.fixture()
def gaussian_full_genecp_outfile(gaussian_outputs_test_directory):
    return os.path.join(
        gaussian_outputs_test_directory,
        "silver_chloride_full_genecp.log",
    )


# Gaussian output file for frozen coordinates
@pytest.fixture()
def gaussian_frozen_opt_outfile(gaussian_outputs_test_directory):
    gaussian_frozen_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "frozen_coordinates_opt.log"
    )
    return gaussian_frozen_opt_outfile


# Gaussian output file for modred
@pytest.fixture()
def gaussian_failed_modred_outfile(gaussian_outputs_test_directory):
    gaussian_modred_outfile = os.path.join(
        gaussian_outputs_test_directory, "cage_free_failed_modred.log"
    )
    return gaussian_modred_outfile


# Gaussian output for scan
@pytest.fixture()
def gaussian_failed_scan_outfile(gaussian_outputs_test_directory):
    gaussian_scan_outfile = os.path.join(
        gaussian_outputs_test_directory, "cationic_failed_scan.log"
    )
    return gaussian_scan_outfile


# Gaussian output file for Hirshfeld charges
@pytest.fixture()
def gaussian_hirshfeld_outfile(gaussian_outputs_test_directory):
    gaussian_hirshfeld_outfile = os.path.join(
        gaussian_outputs_test_directory,
        "oxetane_hirshfeld_sp_smd_n_n-DiMethylFormamide.log",
    )
    return gaussian_hirshfeld_outfile


@pytest.fixture()
def gaussian_rc_hirshfeld_outfile(gaussian_outputs_test_directory):
    gaussian_hirshfeld_outfile = os.path.join(
        gaussian_outputs_test_directory,
        "oxetane_rc_hirshfeld_sp_smd_n_n-DiMethylFormamide.log",
    )
    return gaussian_hirshfeld_outfile


# Gaussian output file with custom (generic) SMD solvent
@pytest.fixture()
def gaussian_smd_generic_outfile(gaussian_outputs_test_directory):
    return os.path.join(
        gaussian_outputs_test_directory,
        "benzoic_acid_opt_sp_smd_generic.log",
    )


@pytest.fixture()
def gaussian_ozone_opt_outfile(gaussian_outputs_test_directory):
    gaussian_ozone_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "ozone.log"
    )
    return gaussian_ozone_opt_outfile


@pytest.fixture()
def gaussian_co2_opt_outfile(gaussian_outputs_test_directory):
    gaussian_co2_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "co2.log"
    )
    return gaussian_co2_opt_outfile


@pytest.fixture()
def gaussian_koh_opt_outfile(gaussian_outputs_test_directory):
    gaussian_koh_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "KOH.log"
    )
    return gaussian_koh_opt_outfile


@pytest.fixture()
def gaussian_koh_linear_opt_outfile(gaussian_outputs_test_directory):
    gaussian_koh_linear_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "KOH_linear.log"
    )
    return gaussian_koh_linear_opt_outfile


@pytest.fixture()
def gaussian_he_opt_outfile(gaussian_outputs_test_directory):
    gaussian_he_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "he.log"
    )
    return gaussian_he_opt_outfile


@pytest.fixture()
def gaussian_acetone_opt_outfile(gaussian_outputs_test_directory):
    gaussian_acetone_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "acetone.log"
    )
    return gaussian_acetone_opt_outfile


@pytest.fixture()
def gaussian_benzene_opt_outfile(gaussian_outputs_test_directory):
    gaussian_benzene_opt_outfile = os.path.join(
        gaussian_outputs_test_directory, "benzene.log"
    )
    return gaussian_benzene_opt_outfile


# Gaussian output file for MP2 calculations
@pytest.fixture()
def gaussian_mp2_outputfile(gaussian_outputs_test_directory):
    gaussian_mp2_outfile = os.path.join(
        gaussian_outputs_test_directory, "water_mp2.log"
    )
    return gaussian_mp2_outfile


# Gaussian output file for (failed) ONIOM calculations
@pytest.fixture()
def gaussian_oniom_outputfile(gaussian_outputs_test_directory):
    gaussian_oniom_outfile = os.path.join(
        gaussian_outputs_test_directory, "failed_oniom_b3lypd3_in_uff.log"
    )
    return gaussian_oniom_outfile


# Gaussian output files for pKa calculations
@pytest.fixture()
def gaussian_pKa_HA_optimization_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_ha_optimization_outputfile = os.path.join(
        gaussian_outputs_test_directory, "5PQ_Me_ts1_no_pd_opt.log"
    )
    return gaussian_pka_ha_optimization_outputfile


@pytest.fixture()
def gaussian_pKa_A_optimization_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_a_optimization_outputfile = os.path.join(
        gaussian_outputs_test_directory, "5PQ_Me_ts1_b_no_pd_opt.log"
    )
    return gaussian_pka_a_optimization_outputfile


@pytest.fixture()
def gaussian_pKa_HA_single_point_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_ha_single_point_outputfile = os.path.join(
        gaussian_outputs_test_directory,
        "5PQ_Me_ts1_no_pd_opt_sp_smd_generic.log",
    )
    return gaussian_pka_ha_single_point_outputfile


@pytest.fixture()
def gaussian_pKa_A_single_point_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_a_single_point_outputfile = os.path.join(
        gaussian_outputs_test_directory,
        "5PQ_Me_ts1_b_no_pd_opt_sp_smd_generic.log",
    )
    return gaussian_pka_a_single_point_outputfile


@pytest.fixture()
def gaussian_pKa_HB_optimization_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_hb_optimization_outputfile = os.path.join(
        gaussian_outputs_test_directory, "collidine-H_opt.log"
    )
    return gaussian_pka_hb_optimization_outputfile


@pytest.fixture()
def gaussian_pKa_B_optimization_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_b_optimization_outputfile = os.path.join(
        gaussian_outputs_test_directory, "collidine_opt.log"
    )
    return gaussian_pka_b_optimization_outputfile


@pytest.fixture()
def gaussian_pKa_HB_single_point_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_hb_single_point_outputfile = os.path.join(
        gaussian_outputs_test_directory, "collidine-H_opt_sp_smd_generic.log"
    )
    return gaussian_pka_hb_single_point_outputfile


@pytest.fixture()
def gaussian_pKa_B_single_point_outputfile(gaussian_outputs_test_directory):
    gaussian_pka_b_single_point_outputfile = os.path.join(
        gaussian_outputs_test_directory, "collidine_opt_sp_smd_generic.log"
    )
    return gaussian_pka_b_single_point_outputfile


# Gaussian pbc input files
@pytest.fixture()
def gaussian_pbc_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "pbc")


@pytest.fixture()
def gaussian_pbc_inputs_test_directory(gaussian_pbc_test_directory):
    return os.path.join(gaussian_pbc_test_directory, "com")


@pytest.fixture()
def gaussian_pbc_1d_inputfile(gaussian_pbc_inputs_test_directory):
    gaussian_pbc_1d_inputfile = os.path.join(
        gaussian_pbc_inputs_test_directory, "neoprene_1d.com"
    )
    return gaussian_pbc_1d_inputfile


# Gaussian PBC output files
@pytest.fixture()
def gaussian_pbc_outputs_test_directory(gaussian_pbc_test_directory):
    return os.path.join(gaussian_pbc_test_directory, "log")


@pytest.fixture()
def gaussian_pbc_2d_outputfile(gaussian_pbc_outputs_test_directory):
    gaussian_pbc_2d_outputfile = os.path.join(
        gaussian_pbc_outputs_test_directory, "graphite_2d_opt.log"
    )
    return gaussian_pbc_2d_outputfile


@pytest.fixture()
def gaussian_pbc_3d_outputfile(gaussian_pbc_outputs_test_directory):
    gaussian_pbc_3d_outputfile = os.path.join(
        gaussian_pbc_outputs_test_directory, "gallium_arsenide_3d.log"
    )
    return gaussian_pbc_3d_outputfile


# text path and associated files
@pytest.fixture()
def txt_path(gaussian_test_directory):
    test_txt_path = os.path.join(gaussian_test_directory, "text")
    return os.path.abspath(test_txt_path)


@pytest.fixture()
def reference_genecp_txt_file_from_api(txt_path):
    return os.path.join(txt_path, "genecp_txt_from_api.txt")


@pytest.fixture()
def genecp_txt_file_from_web(txt_path):
    return os.path.join(txt_path, "test_genecp.txt")


@pytest.fixture()
def gen_txt_file_from_web(txt_path):
    return os.path.join(txt_path, "test_gen.txt")


# Gaussian output file from TDDFT
@pytest.fixture()
def tddft_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "tddft")


@pytest.fixture()
def td_outputfile(tddft_test_directory):
    td_outputfile = os.path.join(
        tddft_test_directory, "tddft_r1s50_gas_radical_anion.log"
    )
    return td_outputfile


# Gaussian cube files
@pytest.fixture()
def cube_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "cubes")


@pytest.fixture()
def spin_cube_file(cube_test_directory):
    spin_cube_file = os.path.join(cube_test_directory, "n2_dens.cube")
    return spin_cube_file


# gaussian yaml files
@pytest.fixture()
def gaussian_yaml_settings_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "project_yaml")


@pytest.fixture()
def gaussian_yaml_settings_gas_solv(gaussian_yaml_settings_directory):
    return os.path.join(gaussian_yaml_settings_directory, "gas_solv.yaml")


@pytest.fixture()
def gaussian_yaml_settings_gas_solv_project_name(
    gaussian_yaml_settings_directory,
):
    return os.path.join(gaussian_yaml_settings_directory, "gas_solv")


@pytest.fixture()
def gaussian_yaml_settings_solv(gaussian_yaml_settings_directory):
    return os.path.join(gaussian_yaml_settings_directory, "solv.yaml")


@pytest.fixture()
def gaussian_yaml_settings_qmmm_project_name(
    gaussian_yaml_settings_directory,
):
    return os.path.join(gaussian_yaml_settings_directory, "qmmm")


# gaussian written files
@pytest.fixture()
def gaussian_written_files_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "written_files")


@pytest.fixture()
def gaussian_written_opt_file(gaussian_written_files_directory):
    return os.path.join(gaussian_written_files_directory, "gaussian_opt.com")


@pytest.fixture()
def gaussian_written_pm6_opt_file(gaussian_written_files_directory):
    return os.path.join(
        gaussian_written_files_directory, "gaussian_pm6_opt.com"
    )


@pytest.fixture()
def gaussian_written_opt_file_with_route(gaussian_written_files_directory):
    return os.path.join(
        gaussian_written_files_directory, "gaussian_opt_with_route.com"
    )


@pytest.fixture()
def gaussian_written_modred_file(gaussian_written_files_directory):
    return os.path.join(
        gaussian_written_files_directory, "gaussian_modred.com"
    )


@pytest.fixture()
def gaussian_written_scan_single_degree_of_freedom_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_scan_single_degree_of_freedom.com",
    )


@pytest.fixture()
def gaussian_written_scan_multiple_degrees_of_freedom_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_scan_multiple_degrees_of_freedom.com",
    )


@pytest.fixture()
def gaussian_written_scan_multiple_degrees_of_freedom_with_constraints_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_scan_multiple_degrees_of_freedom_with_constraints.com",
    )


@pytest.fixture()
def gaussian_written_ts_file(gaussian_written_files_directory):
    return os.path.join(gaussian_written_files_directory, "gaussian_ts.com")


@pytest.fixture()
def gaussian_written_qmmm_file(gaussian_written_files_directory):
    return os.path.join(gaussian_written_files_directory, "gaussian_qmmm.com")


@pytest.fixture()
def gaussian_written_qmmm_log_file(gaussian_written_files_directory):
    return os.path.join(
        gaussian_written_files_directory, "gaussian_qmmm_from_log.com"
    )


@pytest.fixture()
def gaussian_written_ts_from_nhc_singlet_log_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory, "gaussian_ts_from_log.com"
    )


@pytest.fixture()
def gaussian_written_sp_from_nhc_singlet_log_with_solvent_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_sp_from_log_with_solvent.com",
    )


@pytest.fixture()
def gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_sp_from_log_with_custom_solvent.com",
    )


@pytest.fixture()
def gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_sp_from_log_with_custom_basis.com",
    )


@pytest.fixture()
def gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_sp_from_log_with_custom_basis_from_api.com",
    )


@pytest.fixture()
def gaussian_modred_with_custom_basis_for_all_atoms_from_api(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_modred_with_custom_basis_for_all_atoms_from_api.com",
    )


@pytest.fixture()
def gaussian_written_opt_from_graphite_2d_pbc_log(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory, "graphite_2d_opt_from_log.com"
    )


@pytest.fixture()
def qmmm_written_xyz_file(gaussian_written_files_directory):
    return os.path.join(gaussian_written_files_directory, "qmmm_written.xyz")


@pytest.fixture()
def qmmm_written_xyz_only_file(gaussian_written_files_directory):
    return os.path.join(
        gaussian_written_files_directory, "qmmm_written_xyz_only.xyz"
    )


# Gaussian folder for thermochemistry analysis
@pytest.fixture()
def gaussian_thermochem_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "thermochem")


@pytest.fixture()
def gaussian_co2_pressure1p5_outfile(gaussian_thermochem_test_directory):
    gaussian_co2_pressure1p5_outfile = os.path.join(
        gaussian_thermochem_test_directory, "co2_pressure1p5.log"
    )
    return gaussian_co2_pressure1p5_outfile


@pytest.fixture()
def gaussian_co2_pressure3_outfile(gaussian_thermochem_test_directory):
    gaussian_co2_pressure3_outfile = os.path.join(
        gaussian_thermochem_test_directory, "co2_pressure3.log"
    )
    return gaussian_co2_pressure3_outfile


# Gaussian folder for boltzmann weighting
@pytest.fixture()
def gaussian_boltzmann_test_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "boltzmann")


@pytest.fixture()
def gaussian_conformer1_outfile(gaussian_boltzmann_test_directory):
    gaussian_conformer1_outfile = os.path.join(
        gaussian_boltzmann_test_directory, "udc3_mCF3_monomer_c1.log"
    )
    return gaussian_conformer1_outfile


@pytest.fixture()
def gaussian_conformer2_outfile(gaussian_boltzmann_test_directory):
    gaussian_conformer2_outfile = os.path.join(
        gaussian_boltzmann_test_directory, "udc3_mCF3_monomer_c4.log"
    )
    return gaussian_conformer2_outfile


# text path and associated files


@pytest.fixture()
def smd_TBME_solvent_parameters_text_file(txt_path):
    return os.path.join(txt_path, "smd_TBME.txt")


@pytest.fixture()
def Ni_def2tzvp_PCHOSi_svp_text_file(txt_path):
    return os.path.join(txt_path, "Ni_def2tzvp_PCHOSi_svp.txt")


############ Orca Fixtures ##################
# master orca test directory
@pytest.fixture()
def orca_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "ORCATests")


# orca input files path and associated files
@pytest.fixture()
def inpfile_path(orca_test_directory):
    test_inpfile_path = os.path.join(orca_test_directory, "inputs")
    return os.path.abspath(test_inpfile_path)


# specific input files
@pytest.fixture()
def water_sp_input_path(inpfile_path):
    return os.path.join(inpfile_path, "water_sp.inp")


@pytest.fixture()
def water_opt_input_path(inpfile_path):
    return os.path.join(inpfile_path, "water_opt.inp")


# orca input files path and associated files
@pytest.fixture()
def orca_inputs_directory(orca_test_directory):
    orca_inputs_directory = os.path.join(orca_test_directory, "inputs")
    return os.path.abspath(orca_inputs_directory)


@pytest.fixture()
def orca_inputs_xyz_directory(orca_inputs_directory):
    """Returns the absolute path to the
    orca inputs that specifies xyz files."""
    orca_inputs_xyz_directory = os.path.join(orca_inputs_directory, "xyz")
    return os.path.abspath(orca_inputs_xyz_directory)


@pytest.fixture()
def orca_input_nebts_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca
    input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "neb_TS_rot1.inp")


@pytest.fixture()
def orca_input_nebts_reactant_xyz_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca
    input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "R-1a_opt.xyz")


@pytest.fixture()
def orca_input_nebts_product_xyz_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca
    input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "S-1a_opt.xyz")


@pytest.fixture()
def orca_input_nebts_ts_xyz_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca
    input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "TS_rot1.xyz")


@pytest.fixture()
def orca_epr_solv(orca_inputs_directory):
    return os.path.join(orca_inputs_directory, "ORCA_Test_0829.inp")


@pytest.fixture()
def orca_faulty_solv(orca_inputs_directory):
    return os.path.join(orca_inputs_directory, "faulty_solv.inp")


@pytest.fixture()
def orca_outputs_directory(orca_test_directory):
    orca_outputs_directory = os.path.join(orca_test_directory, "outputs")
    return os.path.abspath(orca_outputs_directory)


@pytest.fixture()
def water_sp_gas_path(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "water_dlpno_ccsdt_sp.out")


@pytest.fixture()
def water_output_gas_path(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "water_opt.out")


@pytest.fixture()
def orca_he_output_freq(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "He_freq.out")


@pytest.fixture()
def orca_co2_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "CO2.out")


@pytest.fixture()
def orca_koh_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "KOH.out")


@pytest.fixture()
def orca_sn2_ts_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "sn2_ts.out")


@pytest.fixture()
def dlpno_ccsdt_sp_full_print(orca_outputs_directory):
    return os.path.join(
        orca_outputs_directory, "dlpno_ccsdt_singlepoint_neutral_in_cpcm.out"
    )


@pytest.fixture()
def hirshfeld_full_print(orca_outputs_directory):
    return os.path.join(
        orca_outputs_directory, "udc3_ts1_c15_sp_hirshfeld.out"
    )


@pytest.fixture()
def fe2_singlet_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "fe2_singlet.out")


@pytest.fixture()
def fe2_triplet_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "fe2_triplet.out")


@pytest.fixture()
def fe2_quintet_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "fe2_quintet.out")


@pytest.fixture()
def fe3_doublet_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "fe3_doublet.out")


@pytest.fixture()
def fe3_quartet_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "fe3_quartet.out")


@pytest.fixture()
def fe3_sextet_output(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "fe3_sextet.out")


@pytest.fixture()
def water_engrad_path(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "water_opt.engrad")


@pytest.fixture()
def orca_fixed_atoms(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "phenol_fixed_atoms.out")


@pytest.fixture()
def orca_fixed_bonds_and_angles(orca_outputs_directory):
    return os.path.join(
        orca_outputs_directory, "phenol_fixed_bond_and_angles.out"
    )


@pytest.fixture()
def orca_fixed_dihedral(orca_outputs_directory):
    return os.path.join(
        orca_outputs_directory, "phenylalanine_fixed_dihedral.out"
    )


@pytest.fixture()
def orca_two_layer_qmmmm_output_file(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "methanol_ethane_qmmm.out")


@pytest.fixture()
def orca_neb_output_file(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "neb_R-TS1-Si.out")


@pytest.fixture()
def orca_errors_directory(orca_test_directory):
    orca_errors_directory = os.path.join(orca_test_directory, "error_files")
    return os.path.abspath(orca_errors_directory)


@pytest.fixture()
def gtoint_errfile(orca_errors_directory):
    return os.path.join(orca_errors_directory, "GTOInt_error.out")


# orca written files
@pytest.fixture()
def orca_written_files_directory(orca_test_directory):
    orca_written_files = os.path.join(orca_test_directory, "written_files")
    return orca_written_files


@pytest.fixture()
def orca_written_opt_file(orca_written_files_directory):
    return os.path.join(orca_written_files_directory, "orca_opt.inp")


@pytest.fixture()
def orca_written_opt_file_with_route(orca_written_files_directory):
    return os.path.join(
        orca_written_files_directory, "orca_opt_with_route.inp"
    )


@pytest.fixture()
def orca_written_modred_file(orca_written_files_directory):
    return os.path.join(orca_written_files_directory, "orca_modred.inp")


@pytest.fixture()
def orca_written_scan_single_degree_of_freedom_file(
    orca_written_files_directory,
):
    return os.path.join(
        orca_written_files_directory, "orca_scan_single_degree_of_freedom.inp"
    )


@pytest.fixture()
def orca_written_scan_multiple_degrees_of_freedom_file(
    orca_written_files_directory,
):
    return os.path.join(
        orca_written_files_directory,
        "orca_scan_multiple_degrees_of_freedom.inp",
    )


@pytest.fixture()
def orca_written_scan_multiple_degrees_of_freedom_with_constraints_file(
    orca_written_files_directory,
):
    return os.path.join(
        orca_written_files_directory,
        "orca_scan_multiple_degrees_of_freedom_with_constraints.inp",
    )


@pytest.fixture()
def orca_written_ts_file(orca_written_files_directory):
    return os.path.join(orca_written_files_directory, "orca_ts.inp")


@pytest.fixture()
def orca_written_ts_from_nhc_singlet_log_file(orca_written_files_directory):
    return os.path.join(orca_written_files_directory, "orca_ts_from_log.inp")


@pytest.fixture()
def orca_written_sp_from_nhc_singlet_log_with_solvent_file(
    orca_written_files_directory,
):
    return os.path.join(
        orca_written_files_directory, "orca_sp_from_log_with_solvent.inp"
    )


@pytest.fixture()
def orca_written_he_monoatomic_opt_file(orca_written_files_directory):
    return os.path.join(
        orca_written_files_directory, "orca_he_monoatomic_opt.inp"
    )


@pytest.fixture()
def orca_written_neb_file(orca_written_files_directory):
    return os.path.join(orca_written_files_directory, "orca_neb_TS_rot1.inp")


# orca yaml files
@pytest.fixture()
def orca_yaml_settings_directory(orca_test_directory):
    return os.path.join(orca_test_directory, "project_yaml")


@pytest.fixture()
def orca_yaml_settings_gas_solv_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "gas_solv")


@pytest.fixture()
def orca_yaml_settings_orca_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "orca")


@pytest.fixture()
def orca_yaml_settings_custom_solv_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "custom_solv")


@pytest.fixture()
def orca_yaml_settings_custom_solv_cosmors_project_name(
    orca_yaml_settings_directory,
):
    return os.path.join(orca_yaml_settings_directory, "custom_solv_cosmors")


# master xTB test directory
@pytest.fixture()
def xtb_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "XTBTests")


@pytest.fixture()
def xtb_inputs_directory(xtb_test_directory):
    xtb_inputs_directory = os.path.join(xtb_test_directory, "inputs")
    return os.path.abspath(xtb_inputs_directory)


@pytest.fixture()
def xtb_default_inputfile(xtb_inputs_directory):
    return os.path.join(xtb_inputs_directory, "default.inp")


@pytest.fixture()
def xtb_sp_alpb_inputfile(xtb_inputs_directory):
    return os.path.join(xtb_inputs_directory, "alpb_water.inp")


@pytest.fixture()
def xtb_outputs_directory(xtb_test_directory):
    xtb_outputs_directory = os.path.join(xtb_test_directory, "outputs")
    return os.path.abspath(xtb_outputs_directory)


@pytest.fixture()
def xtb_co2_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "co2_ohess")


@pytest.fixture()
def xtb_water_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "water_ohess")


@pytest.fixture()
def xtb_cyclopentadienyl_anion_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "cyclopentadienyl_anion_opt")


@pytest.fixture()
def xtb_p_benzyne_opt_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "p_benzyne_opt_alpb_toluene")


@pytest.fixture()
def xtb_p_benzyne_sp_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "p_benzyne_sp_alpb_toluene")


@pytest.fixture()
def xtb_acetaldehyde_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "acetaldehyde_hess")


@pytest.fixture()
def xtb_he_outfolder(xtb_outputs_directory):
    return os.path.join(xtb_outputs_directory, "he_hess")


# test for structure.py
@pytest.fixture()
def structure_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "StructuresTests")


@pytest.fixture()
def xyz_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "xyz")


@pytest.fixture()
def single_molecule_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "crest_best.xyz")


@pytest.fixture()
def multiple_molecules_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "crest_conformers.xyz")


@pytest.fixture()
def xtb_optimized_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "ts_xtbopt.xyz")


@pytest.fixture()
def chemsmart_generated_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "frozen_coordinates_opt.xyz")


@pytest.fixture()
def extended_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "crystal.extxyz")


@pytest.fixture()
def chemdraw_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "chemdraw")


@pytest.fixture()
def single_molecule_cdxml_file_benzene(chemdraw_directory):
    return os.path.join(chemdraw_directory, "benzene.cdxml")


@pytest.fixture()
def single_molecule_cdxml_file_methane(chemdraw_directory):
    return os.path.join(chemdraw_directory, "methane.cdxml")


@pytest.fixture()
def multi_molecule_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "two_molecules.cdxml")


@pytest.fixture()
def single_molecule_cdx_file_imidazole(chemdraw_directory):
    return os.path.join(chemdraw_directory, "imidazole.cdx")


@pytest.fixture()
def complex_molecule_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "complex_molecule.cdxml")


@pytest.fixture()
def metal_ligand_molecules_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "metal_ligands.cdxml")


@pytest.fixture()
def colored_proton_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "phenol.cdxml")


@pytest.fixture()
def colored_implicit_proton_cdxml_file(chemdraw_directory):
    """Returns the path to a CDXML file with a colored implicit proton, which should be treated as a colored explicit proton."""
    return os.path.join(chemdraw_directory, "phenol_implicit_proton.cdxml")


@pytest.fixture()
def colored_proton_two_molecule_cdxml_file(chemdraw_directory):
    return os.path.join(chemdraw_directory, "phenol_two_molecule.cdxml")


@pytest.fixture()
def pka_scale_cdxml_file(chemdraw_directory):
    """Multi-fragment CDXML with nested groups and coloured acidic protons."""
    return os.path.join(chemdraw_directory, "pka_scale.cdxml")


@pytest.fixture()
def utils_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "YAMLTests")


@pytest.fixture()
def server_yaml_file(utils_test_directory):
    return os.path.join(utils_test_directory, "server.yaml")


### Server and JobRunner fixtures


@pytest.fixture()
def gaussian_project_config_dir(tmp_path):
    """Minimal Gaussian project config under a temporary CHEMSMART config root."""
    config_root = tmp_path / "chemsmart_cfg"
    gaussian_cfg = config_root / "gaussian"
    gaussian_cfg.mkdir(parents=True)
    (gaussian_cfg / "test.yaml").write_text(
        "gas:\n  functional: B3LYP\n  basis: def2-SVP\n"
        "solv:\n  functional: B3LYP\n  basis: def2-SVP\n"
        "  solvent_model: smd\n  solvent_id: water\n"
    )
    return config_root


@pytest.fixture()
def pbs_server(server_yaml_file):
    return Server.from_yaml(server_yaml_file)


@pytest.fixture()
def gaussian_jobrunner_no_scratch(pbs_server):
    return FakeGaussianJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def gaussian_jobrunner_scratch(tmpdir, pbs_server):
    return FakeGaussianJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def orca_jobrunner_no_scratch(pbs_server):
    return FakeORCAJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def nciplot_jobrunner_no_scratch(pbs_server):
    return FakeNCIPLOTJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def nciplot_jobrunner_scratch(tmpdir, pbs_server):
    return FakeNCIPLOTJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


## pytest fixtures for molecules


@pytest.fixture()
def methyl3hexane_molecule():
    symbols = [
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "C",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
        "H",
    ]
    coords = np.array(
        [
            [0.828, -0.5939, -0.4105],
            [-0.6074, -0.1292, -0.7341],
            [1.6188, 0.3738, 0.4896],
            [-1.5335, 0.0189, 0.4761],
            [0.8379, -1.9975, 0.2043],
            [1.7539, 1.7749, -0.0888],
            [-2.8977, 0.553, 0.0634],
            [1.3601, -0.6621, -1.3694],
            [-1.0561, -0.849, -1.4317],
            [-0.5655, 0.8214, -1.2791],
            [2.6287, -0.0297, 0.6378],
            [1.1655, 0.4353, 1.4856],
            [-1.0999, 0.707, 1.2088],
            [-1.6692, -0.9475, 0.973],
            [1.8635, -2.3744, 0.2836],
            [0.4049, -2.004, 1.2096],
            [0.2737, -2.7005, -0.4174],
            [2.1704, 1.7451, -1.1005],
            [2.4243, 2.3747, 0.5353],
            [0.7895, 2.2899, -0.1235],
            [-3.384, -0.1214, -0.6486],
            [-2.8085, 1.5396, -0.4023],
            [-3.5477, 0.6486, 0.9389],
        ]
    )
    methyl3hexane = Molecule(symbols=symbols, positions=coords)
    return methyl3hexane


@pytest.fixture()
def constrained_atoms():
    """Fixture to create a simple Ar2 dimer with constraints."""
    from ase import Atoms
    from ase.calculators.lj import LennardJones
    from ase.constraints import FixAtoms, FixBondLength

    # Simple Ar2 dimer with a reasonable separation
    r0 = 3.5  # Å
    atoms = Atoms(
        "Ar2", positions=[(0.0, 0.0, 0.0), (r0, 0.0, 0.0)], pbc=False
    )

    # Light-weight calculator for tests
    atoms.calc = LennardJones()  # defaults are fine for unit tests

    # Constraints:
    #  - Fix the first atom in space
    #  - Keep the Ar–Ar bond length fixed at its initial value
    constraints = [
        FixAtoms(indices=[0]),
        FixBondLength(0, 1),
    ]
    # set the constraints on the Atoms object
    atoms.set_constraint(constraint=constraints)

    # set velocity
    atoms.set_velocities([[0, 0, 0], [0, 0, 0]])  # Set zero velocities

    return atoms


@pytest.fixture()
def io_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "IOTests")


@pytest.fixture()
def constrained_pbc_db_file(io_test_directory):
    """Fixture of a .db file containing constrained PBC database
    from heterogeneous catalysis."""
    return os.path.join(
        io_test_directory, "heterogenous_pbc_constraints_5images.db"
    )


## fixtures for mixins


# pytest fixtures for Popen


# Use built-in caplog fixture for capturing log messages
@pytest.fixture()
def capture_log(caplog):
    """
    Fixture to capture log messages.

    Captures messages from the root logger at DEBUG level by default.
    """
    caplog.set_level(logging.DEBUG, logger="")  # "" for root logger
    return caplog


############ Iterate Fixtures ##################


# ── InChIKey test data ──
@pytest.fixture()
def inchikey_test_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "inchikey")


@pytest.fixture()
def inchikey_normal_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "normal_testing",
        "inchikey_normal_testing.xyz",
    )


@pytest.fixture()
def inchikey_r_enantiomer_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "enantiomer_testing",
        "inchikey_r_enantiomer.xyz",
    )


@pytest.fixture()
def inchikey_s_enantiomer_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "enantiomer_testing",
        "inchikey_s_enantiomer.xyz",
    )


@pytest.fixture()
def inchikey_large_molecule_c3_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "large_molecule_testing",
        "inchikey_large_molecule_c3.xyz",
    )


@pytest.fixture()
def inchikey_large_molecule_c2_file(inchikey_test_directory):
    return os.path.join(
        inchikey_test_directory,
        "large_molecule_testing",
        "inchikey_large_molecule_c2.xyz",
    )


# ── CXSMILES test data ──
@pytest.fixture()
def cxsmiles_test_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "cxsmiles")


@pytest.fixture()
def cxsmiles_normal_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "normal_testing",
        "rotamer_normal_testing.xyz",
    )


@pytest.fixture()
def cxsmiles_r_enantiomer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "enantiomer_testing",
        "rotamer_r_enantiomer.xyz",
    )


@pytest.fixture()
def cxsmiles_s_enantiomer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "enantiomer_testing",
        "rotamer_s_enantiomer.xyz",
    )


@pytest.fixture()
def cxsmiles_r_rotamer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "rotamer_testing",
        "cxsmiles_r_rotamer.xyz",
    )


@pytest.fixture()
def cxsmiles_s_rotamer_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "rotamer_testing",
        "cxsmiles_s_rotamer.xyz",
    )


@pytest.fixture()
def cxsmiles_large_molecule_c2_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "rotamer_large_molecule_c2.xyz",
    )


@pytest.fixture()
def cxsmiles_large_molecule_c3_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "rotamer_large_molecule_c3.xyz",
    )


@pytest.fixture()
def cxsmiles_expected_large_c2_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "expected_cxsmiles_c2.txt",
    )


@pytest.fixture()
def cxsmiles_expected_large_c3_file(cxsmiles_test_directory):
    return os.path.join(
        cxsmiles_test_directory,
        "large_molecule_testing",
        "expected_cxsmiles_c3.txt",
    )


############ Molecule Fixtures for RDKit / PDB conversion tests ##################


############ Canonical Geometry / Structure ID Fixtures ##################
@pytest.fixture()
def canonical_test_directory(structure_test_directory):
    return os.path.join(structure_test_directory, "canonical")


@pytest.fixture()
def canonical_formaldehyde_file(canonical_test_directory):
    """Formaldehyde (CH2O) in its reference orientation — planar, C2v symmetry."""
    return os.path.join(canonical_test_directory, "formaldehyde.xyz")


@pytest.fixture()
def canonical_formaldehyde_trans_rot_file(canonical_test_directory):
    """Formaldehyde translated and rotated from canonical_formaldehyde_file."""
    return os.path.join(canonical_test_directory, "formaldehyde_trans_rot.xyz")


@pytest.fixture()
def canonical_formaldehyde_perturbed_file(canonical_test_directory):
    """Formaldehyde with coordinates perturbed at ~1e-7 Å level from canonical_formaldehyde_file."""
    return os.path.join(canonical_test_directory, "formaldehyde_perturbed.xyz")


@pytest.fixture()
def canonical_methane_file(canonical_test_directory):
    """Methane (CH4) in its reference orientation — tetrahedral, Td symmetry."""
    return os.path.join(canonical_test_directory, "methane.xyz")


@pytest.fixture()
def canonical_methane_trans_rot_file(canonical_test_directory):
    """Methane translated and rotated from canonical_methane_file."""
    return os.path.join(canonical_test_directory, "methane_trans_rot.xyz")


@pytest.fixture()
def canonical_methane_distorted_file(canonical_test_directory):
    """Methane with one C-H bond elongated by ~2e-3 Å from canonical_methane_file."""
    return os.path.join(canonical_test_directory, "methane_distorted.xyz")


@pytest.fixture()
def canonical_3b_file(canonical_test_directory):
    """3b (C17H17NOS) in its reference orientation — large, low-symmetry (C1) molecule."""
    return os.path.join(canonical_test_directory, "3b.xyz")


@pytest.fixture()
def canonical_3b_trans_rot_file(canonical_test_directory):
    """3b translated and rotated from canonical_3b_file."""
    return os.path.join(canonical_test_directory, "3b_trans_rot.xyz")


@pytest.fixture()
def canonical_3b_permuted_file(canonical_test_directory):
    """3b with input atom order permuted — same coordinates as canonical_3b_file."""
    return os.path.join(canonical_test_directory, "3b_permuted.xyz")


@pytest.fixture()
def canonical_r_bromochlorofluoromethane_file(canonical_test_directory):
    """(R)-Bromochlorofluoromethane — one enantiomer of a chiral CHBrClF molecule."""
    return os.path.join(
        canonical_test_directory, "R-Bromochlorofluoromethane.xyz"
    )


@pytest.fixture()
def canonical_s_bromochlorofluoromethane_file(canonical_test_directory):
    """(S)-Bromochlorofluoromethane — non-superimposable mirror image of the R enantiomer."""
    return os.path.join(
        canonical_test_directory, "S-Bromochlorofluoromethane.xyz"
    )


############ Database Fixtures ##################
