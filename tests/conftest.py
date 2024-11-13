import os
import pytest


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
    return os.path.abspath("data")


# master gaussian test directory
@pytest.fixture()
def gaussian_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "GaussianTests")


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


############ Orca Fixtures ##################
# master orca test directory
@pytest.fixture()
def orca_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "ORCATests")


# orca input files path and associated files
@pytest.fixture()
def inpfile_path(orca_test_directory):
    test_inpfile_path = os.path.join(orca_test_directory, "orca_inputs")
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
    orca_inputs_directory = os.path.join(orca_test_directory, "orca_inputs")
    return os.path.abspath(orca_inputs_directory)


@pytest.fixture()
def orca_dias_directory(orca_test_directory):
    orca_dias_directory = os.path.join(orca_test_directory, "orca_dias")
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
    orca_outputs_directory = os.path.join(orca_test_directory, "orca_outputs")
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
def water_engrad_path(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "water_opt.engrad")


@pytest.fixture()
def orca_errors_directory(orca_test_directory):
    orca_errors_directory = os.path.join(
        orca_test_directory, "orca_error_files"
    )
    return os.path.abspath(orca_errors_directory)


@pytest.fixture()
def gtoint_errfile(orca_errors_directory):
    return os.path.join(orca_errors_directory, "GTOInt_error.out")


# orca yaml files
@pytest.fixture()
def orca_yaml_settings_directory(orca_test_directory):
    return os.path.join(orca_test_directory, "yaml_settings")


@pytest.fixture()
def orca_yaml_settings_defaults(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "defaults.yaml")


@pytest.fixture()
def orca_yaml_settings_gas_solv(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "gas_solv.yaml")


@pytest.fixture()
def orca_yaml_settings_solv(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "solv.yaml")
