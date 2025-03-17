import os

import pytest
import rdkit.Chem.rdDistGeom as rdDistGeom
from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.runner import FakeGaussianJobRunner
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
def gaussian_written_scan_file(gaussian_written_files_directory):
    return os.path.join(gaussian_written_files_directory, "gaussian_scan.com")


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
    orca_errors_directory = os.path.join(orca_test_directory, "error_files")
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
def jobrunner_no_scratch(pbs_server):
    return FakeGaussianJobRunner(server=pbs_server, scratch=False, fake=True)


@pytest.fixture()
def jobrunner_scratch(pbs_server):
    return FakeGaussianJobRunner(server=pbs_server, scratch=True, fake=True)


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


## xtb fixtures
# master xtb test directory
@pytest.fixture()
def xtb_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "XTBTests")


@pytest.fixture()
def xtb_inputs_directory(xtb_test_directory):
    return os.path.join(xtb_test_directory, "inputs")


@pytest.fixture()
def xtb_default_inputfile(xtb_inputs_directory):
    xtb_default_inputfile = os.path.join(xtb_inputs_directory, "default.inp")
    return xtb_default_inputfile


@pytest.fixture()
def xtb_outputs_directory(xtb_test_directory):
    return os.path.join(xtb_test_directory, "outputs")


@pytest.fixture()
def xtb_sp_outfile(xtb_outputs_directory):
    xtb_sp_outfile = os.path.join(xtb_outputs_directory, "water_sp.out")
    return xtb_sp_outfile


@pytest.fixture()
def xtb_opt_outfile(xtb_outputs_directory):
    xtb_opt_outfile = os.path.join(xtb_outputs_directory, "water_opt.out")
    return xtb_opt_outfile


@pytest.fixture()
def xtb_opt_gbsa_outfile(xtb_outputs_directory):
    xtb_opt_gbsa_outfile = os.path.join(
        xtb_outputs_directory, "pyridine_opt_acetonitrile.out"
    )
    return xtb_opt_gbsa_outfile


@pytest.fixture()
def xtb_hess_outfile(xtb_outputs_directory):
    xtb_hess_outfile = os.path.join(
        xtb_outputs_directory, "pyridine_hess_acetonitrile.out"
    )
    return xtb_hess_outfile
