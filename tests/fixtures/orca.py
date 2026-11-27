"""ORCA pKa settings and data-file fixtures."""

import os

import pytest


############ ORCA pKa Settings Fixtures ##################
@pytest.fixture()
def orca_reference_xyz_file(tmp_path, ethanol_molecule):
    """Ethanol written to ``tmp_path/href.xyz`` as a pKa reference file."""
    path = tmp_path / "href.xyz"
    ethanol_molecule.write(str(path))
    return str(path)


@pytest.fixture()
def orca_pka_settings():
    """Real ``ORCApKaJobSettings`` (ethanol O-H proton) without a reference."""
    from chemsmart.jobs.orca.settings import ORCApKaJobSettings

    return ORCApKaJobSettings(proton_index=4, charge=0, multiplicity=1)


@pytest.fixture()
def orca_pka_settings_with_reference(orca_reference_xyz_file):
    """``ORCApKaJobSettings`` with an on-disk ethanol reference file."""
    from chemsmart.jobs.orca.settings import ORCApKaJobSettings

    return ORCApKaJobSettings(
        proton_index=4,
        charge=0,
        multiplicity=1,
        reference_file=orca_reference_xyz_file,
        reference_proton_index=4,
        reference_charge=0,
        reference_multiplicity=1,
    )


@pytest.fixture()
def orca_td_settings(orca_td_project_name):
    """``ORCATDDFTJobSettings`` built from the ``orca_td`` project's ``td:``
    block, as a closed-shell neutral (charge 0, multiplicity 1) system.

    Each test gets a fresh object, so tests may mutate it freely.
    """
    from chemsmart.jobs.orca.settings import ORCATDDFTJobSettings
    from chemsmart.settings.orca import ORCAProjectSettings

    td_settings = ORCAProjectSettings.from_project(
        orca_td_project_name
    ).td_settings()
    td_settings.charge = 0
    td_settings.multiplicity = 1
    return ORCATDDFTJobSettings(**td_settings.__dict__)


@pytest.fixture()
def write_orca_td_input(
    tmpdir, single_molecule_xyz_file, orca_jobrunner_no_scratch
):
    """Return a callable that writes a real ORCA TD input and returns its text.

    Builds an ``ORCATDDFTJob`` from ``single_molecule_xyz_file`` with the
    given settings, runs the real ``ORCAInputWriter`` into ``tmpdir`` and
    reads back ``<label>.inp``, so tests assert on genuine writer output.

    Usage in tests::

        def test_something(orca_td_settings, write_orca_td_input):
            content = write_orca_td_input(orca_td_settings)
    """
    from chemsmart.jobs.orca import ORCATDDFTJob
    from chemsmart.jobs.orca.writer import ORCAInputWriter

    def _write(settings, label="orca_td"):
        job = ORCATDDFTJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label=label,
            jobrunner=orca_jobrunner_no_scratch,
        )
        ORCAInputWriter(job=job).write(target_directory=str(tmpdir))
        with open(os.path.join(str(tmpdir), f"{label}.inp")) as f:
            return f.read()

    return _write


@pytest.fixture()
def build_orca_td_args(orca_td_project_name):
    """Build ``orca ... td ...`` CLI argv, defaulting to the ``td`` project.

    ``pre_td`` holds orca-group options placed before ``td`` (e.g.
    ``-r Opt`` or ``-i 1-2``); ``project_name`` overrides the project YAML.

    Usage in tests::

        def test_something(build_orca_td_args, single_molecule_xyz_file):
            args = build_orca_td_args(
                single_molecule_xyz_file, "--nroots", "5",
                pre_td=("-i", "1-2"),
            )
    """

    def _build(
        xyz_file,
        *td_args,
        pre_td=(),
        charge=0,
        multiplicity=1,
        project_name=None,
    ):
        project_name = project_name or orca_td_project_name
        return [
            "-p",
            os.path.basename(project_name),
            "-f",
            xyz_file,
            "-c",
            str(charge),
            "-m",
            str(multiplicity),
            *pre_td,
            "td",
            *td_args,
        ]

    return _build


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
def water_sdf_file(structure_test_directory):
    """Single-molecule SDF (V2000 molfile) for water."""
    return os.path.join(structure_test_directory, "sdf", "water.sdf")


@pytest.fixture()
def water_sdf_file_content(water_sdf_file):
    """Returns the content of the water.sdf file as a string."""
    with open(water_sdf_file, "r") as f:
        content = f.read()
    return content


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
def orca_qmmm_input_file(orca_inputs_directory):
    return os.path.join(orca_inputs_directory, "dna_qmmm.inp")


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
def orca_he_output_freq_new(orca_outputs_directory):
    return os.path.join(orca_outputs_directory, "He_freq_new.out")


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


@pytest.fixture()
def orca_yaml_settings_neb_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "neb")


@pytest.fixture()
def orca_yaml_settings_custom_solv_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "custom_solv")


@pytest.fixture()
def orca_yaml_settings_custom_solv_cosmors_project_name(
    orca_yaml_settings_directory,
):
    return os.path.join(orca_yaml_settings_directory, "custom_solv_cosmors")


@pytest.fixture
def orca_td_project_name(orca_yaml_settings_directory):
    return os.path.join(orca_yaml_settings_directory, "orca_td")


@pytest.fixture
def orca_td_extras_project_name(orca_yaml_settings_directory):
    """Project name whose ``td:`` block sets ``additional_route_parameters``."""
    return os.path.join(orca_yaml_settings_directory, "orca_td_extras")


@pytest.fixture
def orca_td_numfreq_project_name(orca_yaml_settings_directory):
    """Project name whose ``td:`` block requests ``NumFreq`` via YAML."""
    return os.path.join(orca_yaml_settings_directory, "orca_td_numfreq")


@pytest.fixture
def orca_td_inherited_project_name(orca_yaml_settings_directory):
    """Project name whose ``td:`` block omits ``freq``; combined with the
    ``ORCAJobSettings`` defaults this yields ``freq=True`` inherited into the
    TD settings, without any ``additional_route_parameters``."""
    return os.path.join(orca_yaml_settings_directory, "orca_td_inherited")
