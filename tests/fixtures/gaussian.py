"""Gaussian writer, pKa/UV-Vis job and data-file fixtures."""

import os
from unittest.mock import patch

import pytest

from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings
from tests.helpers import (
    ConcreteGaussianProjectSettings,
    FakePKaSettings,
    make_pka_sub_job_factory,
)

############ Gaussian Writer Fixtures ##################


@pytest.fixture()
def make_nciplot_writer(tmp_path):
    """Factory for an ``NCIPLOTInputWriter`` built without ``__init__``.

    Bypasses ``NCIPLOTJob``/``JobRunner`` entirely, since most of the
    writer's validation branches are pure functions of
    ``self.job.filenames``/``self.settings`` and need no real job or
    running directory. The fake job's ``folder`` and the jobrunner's
    ``running_directory`` both point at the test's ``tmp_path``.
    Keyword arguments override attributes on the default
    ``NCIPLOTJobSettings``.
    """
    from types import SimpleNamespace

    from chemsmart.jobs.nciplot.writer import NCIPLOTInputWriter

    def _make(**setting_overrides):
        writer = NCIPLOTInputWriter.__new__(NCIPLOTInputWriter)
        writer.job = SimpleNamespace(
            filenames=None, label="test", folder=str(tmp_path)
        )
        writer.settings = NCIPLOTJobSettings()
        writer.jobrunner = SimpleNamespace(running_directory=str(tmp_path))
        for key, value in setting_overrides.items():
            setattr(writer.settings, key, value)
        return writer

    return _make


############ Gaussian UV-Vis / pKa Job Fixtures ##################


@pytest.fixture
def _clear_shared_reference_cache():
    """Reset ``GaussianpKaJob``'s class-level reference-molecule cache.

    ``_shared_reference_molecule_cache`` is shared by every
    ``GaussianpKaJob`` instance; clearing it around each test stops
    state leaking between tests.
    """
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob

    GaussianpKaJob._shared_reference_molecule_cache.clear()
    yield
    GaussianpKaJob._shared_reference_molecule_cache.clear()


@pytest.fixture()
def pka_settings():
    return FakePKaSettings()


@pytest.fixture()
def pka_settings_with_reference():
    return FakePKaSettings(
        reference_file="href.xyz",
        reference_proton_index=1,
        reference_charge=0,
        reference_multiplicity=1,
    )


@pytest.fixture()
def opt_factory():
    """Recording side_effect factory for pKa's mocked ``GaussianOptJob``."""
    return make_pka_sub_job_factory()


@pytest.fixture()
def sp_factory():
    """Recording side_effect factory for pKa's mocked SP job class."""
    return make_pka_sub_job_factory()


@pytest.fixture()
def patched_sub_jobs(opt_factory, sp_factory):
    """Patch pKa's opt/SP sub-job classes with the recording factories."""
    with (
        patch(
            "chemsmart.jobs.gaussian.pka.GaussianOptJob",
            side_effect=opt_factory,
        ),
        patch(
            "chemsmart.jobs.gaussian.pka.GaussianSinglePointJob",
            side_effect=sp_factory,
        ),
    ):
        yield opt_factory, sp_factory


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
def gaussian_link_irc_outputfile(gaussian_link_outputs_test_directory):
    return os.path.join(
        gaussian_link_outputs_test_directory, "link_irc_output_irc_link.log"
    )


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
def gaussian_dna_link_sp_outputfile(gaussian_link_outputs_test_directory):
    gaussian_link_outfile = os.path.join(
        gaussian_link_outputs_test_directory, "dna_link_sp.log"
    )
    return gaussian_link_outfile


@pytest.fixture()
def gaussian_dppeFeCl2_link_opt_outputfile(
    gaussian_link_outputs_test_directory,
):
    gaussian_link_opt_outfile = os.path.join(
        gaussian_link_outputs_test_directory,
        "dppeFeCl2_opt_quintet_link_opt_link.log",
    )
    return gaussian_link_opt_outfile


@pytest.fixture()
def gaussian_dppeFeCl2_link_opt_failed_outputfile(
    gaussian_link_outputs_test_directory,
):
    gaussian_link_failed_outfile = os.path.join(
        gaussian_link_outputs_test_directory,
        "dppeFeCl2_phenyldioxazolone_opt_triplet_opt_error_termination_link.log",
    )
    return gaussian_link_failed_outfile


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


@pytest.fixture()
def gaussian_link_sp_outfile(gaussian_outputs_test_directory):
    gaussian_link_outfile = os.path.join(
        gaussian_outputs_test_directory, "dna_link_sp.log"
    )
    return gaussian_link_outfile


@pytest.fixture()
def gaussian_link_opt_outfile(gaussian_outputs_test_directory):
    gaussian_link_opt_outfile = os.path.join(
        gaussian_outputs_test_directory,
        "dppeFeCl2_opt_quintet_link_opt_link.log",
    )
    return gaussian_link_opt_outfile


@pytest.fixture()
def gaussian_link_failed_outfile(gaussian_outputs_test_directory):
    gaussian_link_failed_outfile = os.path.join(
        gaussian_outputs_test_directory,
        "dppeFeCl2_phenyldioxazolone_opt_triplet_opt_error_termination_link.log",
    )
    return gaussian_link_failed_outfile


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


@pytest.fixture()
def smd_TBME_solvent_parameters_txt_file(txt_path):
    return os.path.join(txt_path, "smd_TBME.txt")


@pytest.fixture()
def Ni_def2tzvp_PCHOSi_svp_txt_file(txt_path):
    return os.path.join(txt_path, "Ni_def2tzvp_PCHOSi_svp.txt")


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


@pytest.fixture()
def esp_cube_file(cube_test_directory):
    esp_cube_file = os.path.join(cube_test_directory, "n2_esp.cube")
    return esp_cube_file


# gaussian yaml files
@pytest.fixture()
def gaussian_yaml_settings_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "project_yaml")


@pytest.fixture()
def gaussian_yaml_settings_defaults(gaussian_yaml_settings_directory):
    return os.path.join(gaussian_yaml_settings_directory, "defaults.yaml")


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
def gaussian_yaml_settings_qmmm(gaussian_yaml_settings_directory):
    return os.path.join(gaussian_yaml_settings_directory, "qmmm.yaml")


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
def gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file_v2(
    gaussian_written_files_directory,
):
    return os.path.join(
        gaussian_written_files_directory,
        "gaussian_sp_from_log_with_custom_basis_from_api_v2.com",
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
def text_directory(gaussian_test_directory):
    return os.path.join(gaussian_test_directory, "text")


@pytest.fixture()
def genecp_text_file_from_web(text_directory):
    return os.path.join(text_directory, "test_genecp.txt")


@pytest.fixture()
def gen_text_file_from_web(text_directory):
    return os.path.join(txt_path, "test_gen.txt")


@pytest.fixture()
def smd_TBME_solvent_parameters_text_file(txt_path):
    return os.path.join(txt_path, "smd_TBME.txt")


@pytest.fixture()
def Ni_def2tzvp_PCHOSi_svp_text_file(txt_path):
    return os.path.join(txt_path, "Ni_def2tzvp_PCHOSi_svp.txt")


@pytest.fixture()
def concrete_gaussian_project_settings():
    """A fresh hardcoded (non-YAML) :class:`ConcreteGaussianProjectSettings`
    instance, exercising ``GaussianProjectSettings``' default methods."""
    return ConcreteGaussianProjectSettings()
