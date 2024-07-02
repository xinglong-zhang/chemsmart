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
