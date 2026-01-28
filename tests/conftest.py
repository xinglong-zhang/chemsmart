import logging
import os
import tempfile

import pytest
import rdkit.Chem.rdDistGeom as rdDistGeom
import yaml
from pytest_mock import MockerFixture
from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.runner import FakeGaussianJobRunner
from chemsmart.jobs.mol.runner import (
    PyMOLAlignJobRunner,
    PyMOLHybridVisualizationJobRunner,
    PyMOLIRCMovieJobRunner,
    PyMOLMOJobRunner,
    PyMOLMovieJobRunner,
    PyMOLVisualizationJobRunner,
)
from chemsmart.jobs.nciplot.runner import FakeNCIPLOTJobRunner
from chemsmart.jobs.orca.runner import FakeORCAJobRunner
from chemsmart.settings.server import Server

# each test runs on cwd to its temp dir
# @pytest.fixture(autouse=True)
# def go_to_tmpdir(request):
#     # Get the fixture dynamically by its name.
#     tmpdir = request.getfixturevalue("tmpdir")
#     # ensure local test created packages can be imported
#     sys.path.insert(0, str(tmpdir))
#     # Chdir only for the duration of the test.
#     with tmpdir.as_cwd():
#         yield


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


@pytest.fixture()
def sdf_file(test_data_directory):
    return os.path.join(
        test_data_directory, "AtomsWrapperTest", "structure.sdf"
    )


# orca input files path and associated files
@pytest.fixture()
def orca_inputs_directory(orca_test_directory):
    orca_inputs_directory = os.path.join(orca_test_directory, "inputs")
    return os.path.abspath(orca_inputs_directory)


@pytest.fixture()
def orca_inputs_xyz_directory(orca_inputs_directory):
    """Returns the absolute path to the orca inputs that specifies xyz files."""
    orca_inputs_xyz_directory = os.path.join(orca_inputs_directory, "xyz")
    return os.path.abspath(orca_inputs_xyz_directory)


@pytest.fixture()
def orca_input_nebts_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "neb_TS_rot1.inp")


@pytest.fixture()
def orca_input_nebts_reactant_xyz_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "R-1a_opt.xyz")


@pytest.fixture()
def orca_input_nebts_product_xyz_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "S-1a_opt.xyz")


@pytest.fixture()
def orca_input_nebts_ts_xyz_file(orca_inputs_xyz_directory):
    """Returns the absolute path to the orca input file for NEB with TS optimization."""
    return os.path.join(orca_inputs_xyz_directory, "TS_rot1.xyz")


@pytest.fixture()
def orca_dias_directory(orca_test_directory):
    orca_dias_directory = os.path.join(orca_test_directory, "dias")
    return os.path.abspath(orca_dias_directory)


@pytest.fixture()
def water_sp_gas_input_path(orca_inputs_directory):
    return os.path.join(orca_inputs_directory, "water_dlpno_ccsdt_sp.inp")


@pytest.fixture()
def water_sp_solv_input_path(orca_inputs_directory):
    return os.path.join(orca_inputs_directory, "water_dlpno_ccsdt_sp_solv.inp")


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
def water_sp_solv_path(orca_outputs_directory):
    return os.path.join(
        orca_outputs_directory, "water_dlpno_ccsdt_sp_solv.out"
    )


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


# orca yaml files
@pytest.fixture()
def orca_yaml_settings_directory(orca_test_directory):
    return os.path.join(orca_test_directory, "project_yaml")


@pytest.fixture()
def orca_yaml_settings_defaults(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "defaults.yaml")


@pytest.fixture()
def orca_yaml_settings_gas_solv(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "gas_solv.yaml")


@pytest.fixture()
def orca_yaml_settings_solv(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "solv.yaml")


@pytest.fixture()
def orca_yaml_settings_gas_solv_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "gas_solv")


@pytest.fixture()
def orca_yaml_settings_solv_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "solv")


@pytest.fixture()
def orca_yaml_settings_orca_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "orca")


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
def dna_hybrid_visualized_xyz_file(xyz_directory):
    return os.path.join(xyz_directory, "dna_hybrid.xyz")


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
def utils_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "YAMLTests")


@pytest.fixture()
def server_yaml_file(utils_test_directory):
    return os.path.join(utils_test_directory, "server.yaml")


### Server and JobRunner fixtures


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
def orca_jobrunner_scratch(tmpdir, pbs_server):
    return FakeORCAJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def pymol_visualization_jobrunner(pbs_server):
    return PyMOLVisualizationJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_hybrid_visualization_jobrunner(pbs_server):
    return PyMOLHybridVisualizationJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_movie_jobrunner(pbs_server):
    return PyMOLMovieJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def nciplot_jobrunner_no_scratch(pbs_server):
    return FakeNCIPLOTJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def nciplot_jobrunner_scratch(tmpdir, pbs_server):
    return FakeNCIPLOTJobRunner(
        scratch_dir=tmpdir, server=pbs_server, scratch=True, fake=True
    )


@pytest.fixture()
def pymol_align_jobrunner(pbs_server):
    return PyMOLAlignJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_ircmovie_jobrunner(pbs_server):
    return PyMOLIRCMovieJobRunner(server=pbs_server, scratch=False)


@pytest.fixture()
def pymol_mo_jobrunner(pbs_server):
    return PyMOLMOJobRunner(server=pbs_server, scratch=False)


## conformers for testing
@pytest.fixture()
def methanol_molecules():
    # molecules for testing
    # methanol
    methanol = Molecule.from_pubchem(identifier="CO")
    # f = open("methanol.xyz", "w")
    # methanol.write_coordinates(f)

    # rotated methanol
    ase_atoms = methanol.to_ase()
    ase_atoms.rotate(90, [0, 0, 1])
    methanol_rot1 = Molecule.from_ase_atoms(ase_atoms)

    ase_atoms = methanol.to_ase()
    ase_atoms.rotate(20, [1, 1, 1])
    methanol_rot2 = Molecule.from_ase_atoms(ase_atoms)

    methanol_molecules = [methanol, methanol_rot1, methanol_rot2]

    return methanol_molecules


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
def methanol_and_ethanol():
    # molecules for testing
    # methanol
    methanol = Molecule.from_pubchem(identifier="CO")

    # ethanol
    ethanol = Molecule.from_pubchem(identifier="CCO")

    methanol_and_ethanol = [methanol, ethanol]
    return methanol_and_ethanol


@pytest.fixture()
def conformers_from_rdkit():
    """Generate multiple conformers for a complex molecule using RDKit."""
    smiles = "O=C([O-])CCn1c(=O)c(=O)[nH]c2cc([N+](=O)[O-])c(-n3ccc(C=NOCc4ccccc4)c3)cc21"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Generate 3D conformers
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = 0xD06F00D
    ps.numThreads = 10
    conf_ids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=300, params=ps)

    # Ensure each conformer is extracted into its own unique RDKit Mol object
    conformers = []
    for conf_id in conf_ids:
        single_conf_mol = Chem.Mol(mol)  # Copy molecule structure
        single_conf_mol.RemoveAllConformers()  # Remove all existing conformers
        single_conf_mol.AddConformer(
            mol.GetConformer(conf_id), assignId=True
        )  # Add only this conformer
        conformers.append(single_conf_mol)

    # Verify that each conformer contains exactly one conformer
    for conf in conformers:
        assert (
            conf.GetNumConformers() == 1
        ), "Each conformer should contain exactly one conformer."

    # Convert to Molecule instances
    conformers_from_rdkit = [
        Molecule.from_rdkit_mol(conf) for conf in conformers
    ]

    return conformers_from_rdkit


@pytest.fixture()
def io_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "IOTests")


@pytest.fixture()
def excel_file(io_test_directory):
    return os.path.join(io_test_directory, "test.xlsx")


@pytest.fixture()
def constrained_pbc_db_file(io_test_directory):
    """Fixture of a .db file containing constrained PBC database
    from heterogeneous catalysis."""
    return os.path.join(
        io_test_directory, "heterogenous_pbc_constraints_5images.db"
    )


## fixtures for mixins
@pytest.fixture()
def temp_text_file():
    with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
        tmp.write("Line1\nLine2\n")
        tmp_name = tmp.name
    yield tmp_name
    os.remove(tmp_name)


@pytest.fixture()
def dummy_yaml_file():
    class DummyYAMLFile:
        def __init__(self):
            self.filename = "dummy.yaml"
            self.content_lines_string = yaml.dump(
                {"key1": "value1", "key2": "value2"}
            )

        @property
        def yaml_contents_dict(self):
            return yaml.safe_load(self.content_lines_string)

        @property
        def yaml_contents_keys(self):
            return self.yaml_contents_dict.keys()

        @property
        def yaml_contents_values(self):
            return self.yaml_contents_dict.values()

        def yaml_contents_by_key(self, key):
            return self.yaml_contents_dict.get(key)

    return DummyYAMLFile()


@pytest.fixture()
def temp_folder_with_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        file1 = os.path.join(tmpdir, "test1.txt")
        file2 = os.path.join(tmpdir, "test2.log")
        with open(file1, "w") as f:
            f.write("Test file 1")
        with open(file2, "w") as f:
            f.write("Test file 2")
        yield tmpdir, file1, file2


# pytest fixtures for Popen
@pytest.fixture()
def mock_popen(mocker):
    """Fixture to mock subprocess.Popen."""
    return mocker.patch("subprocess.Popen")


@pytest.fixture(scope="session")
def session_mocker(pytestconfig):
    """Session-scoped mocker fixture for patching during the test session."""
    from unittest.mock import MagicMock

    mocker = MockerFixture(pytestconfig)
    mock = MagicMock()
    mocker.patch = mock.patch
    mocker.patch.object = mock.patch.object
    yield mocker
    mocker.resetall()


@pytest.fixture(scope="session")
def tests_logger():
    """Fixture to configure the root logger for tests."""
    logger = logging.getLogger()  # Root logger
    logger.setLevel(logging.INFO)
    logger.handlers = []  # Clear handlers to avoid conflicts
    logger.propagate = True
    # Set environment variable to signal test mode
    os.environ["TEST_MODE"] = "1"
    yield logger
    # Clean up
    logger.handlers = []
    os.environ.pop("TEST_MODE", None)


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
@pytest.fixture()
def iterate_test_directory(test_data_directory):
    """Returns the absolute path to tests/data/IterateTests."""
    return os.path.join(test_data_directory, "IterateTests")


@pytest.fixture()
def iterate_input_directory(iterate_test_directory):
    """Returns the absolute path to tests/data/IterateTests/input."""
    return os.path.join(iterate_test_directory, "input")


@pytest.fixture()
def iterate_expected_output_directory(iterate_test_directory):
    """Returns the absolute path to tests/data/IterateTests/expected_output."""
    return os.path.join(iterate_test_directory, "expected_output")


@pytest.fixture()
def iterate_configs_directory(iterate_test_directory):
    """Returns the absolute path to tests/data/IterateTests/configs."""
    return os.path.join(iterate_test_directory, "configs")


@pytest.fixture()
def iterate_regression_config_file(iterate_configs_directory):
    """Returns the absolute path to tests/data/IterateTests/configs/regression_iterate.cfg."""
    return os.path.join(iterate_configs_directory, "regression_iterate.cfg")


@pytest.fixture()
def iterate_timeout_config_file(iterate_configs_directory):
    """Returns the absolute path to tests/data/IterateTests/configs/timeout_iterate.cfg."""
    return os.path.join(iterate_configs_directory, "timeout_iterate.cfg")


@pytest.fixture()
def iterate_expected_output_file(iterate_expected_output_directory):
    """Returns the absolute path to tests/data/IterateTests/expected_output/regression_iterate_SLSQP_lagrange_multipliers_96_6.xyz."""
    return os.path.join(
        iterate_expected_output_directory,
        "regression_iterate_SLSQP_lagrange_multipliers_96_6.xyz",
    )
