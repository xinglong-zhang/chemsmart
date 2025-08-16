import os
from shutil import copy

import numpy as np

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian import (
    GaussianModredJob,
    GaussianOptJob,
    GaussianScanJob,
    GaussianSinglePointJob,
    GaussianTSJob,
)
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.settings.gaussian import GaussianProjectSettings


class TestGaussianJobs:
    def test_gaussian_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file,
    ):
        # get project settings
        # from tests/data/GaussianTests/project_yaml/gas_solv.yaml
        # functional: "b3lyp empiricaldispersion=gd3bj"
        # basis: def2svp
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianOptJob)
        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)

        outputfile = os.path.join(tmpdir, "gaussian_opt_fake.log")

        g16_fake_opt_out = Gaussian16Output(outputfile)
        assert (
            "b3lyp empiricaldispersion=gd3bj def2svp"
            in g16_fake_opt_out.route_string
        )
        mol = Molecule.from_filepath(outputfile)
        assert np.allclose(
            mol.positions[0], [-1.04402, -2.39212, -1.17658], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
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
            "C",
            "O",
            "O",
            "H",
            "C",
            "N",
            "N",
            "N",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "C",
            "C",
            "C",
            "Cl",
            "C",
            "Cl",
            "H",
            "H",
            "O",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "C",
            "H",
            "Cl",
            "C",
            "H",
            "H",
            "H",
        ]

    def test_gaussian_semiempirical_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_pm6_opt_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.semiempirical = "PM6"
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_pm6_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianOptJob)
        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        outputfile = os.path.join(tmpdir, "gaussian_pm6_opt_fake.log")

        g16_fake_opt_out = Gaussian16Output(outputfile)
        assert "pm6" in g16_fake_opt_out.route_string
        mol = Molecule.from_filepath(outputfile)
        assert np.allclose(
            mol.positions[0], [-1.04402, -2.39212, -1.17658], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
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
            "C",
            "O",
            "O",
            "H",
            "C",
            "N",
            "N",
            "N",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "C",
            "C",
            "C",
            "Cl",
            "C",
            "Cl",
            "H",
            "H",
            "O",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "C",
            "H",
            "Cl",
            "C",
            "H",
            "H",
            "H",
        ]

    def test_gaussian_opt_job_with_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file_with_route,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.route_to_be_written = "#p pbepbe/6-31g(d) opt"
        settings.title = "Optimisation job with supplied route"
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianOptJob)
        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        outputfile = os.path.join(tmpdir, "gaussian_opt_fake.log")

        g16_fake_opt_out = Gaussian16Output(outputfile)
        assert "#p pbepbe/6-31g(d) opt" in g16_fake_opt_out.route_string
        mol = Molecule.from_filepath(outputfile)
        assert np.allclose(
            mol.positions[0], [-1.04402, -2.39212, -1.17658], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
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
            "C",
            "O",
            "O",
            "H",
            "C",
            "N",
            "N",
            "N",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "C",
            "C",
            "C",
            "Cl",
            "C",
            "Cl",
            "H",
            "H",
            "O",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "C",
            "H",
            "Cl",
            "C",
            "H",
            "H",
            "H",
        ]

    def test_gaussian_modred_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_modred_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.modred_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = [[1, 2], [3, 4, 5]]
        job = GaussianModredJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_modred",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianModredJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        outputfile = os.path.join(tmpdir, "gaussian_modred_fake.log")

        g16_fake_opt_out = Gaussian16Output(outputfile)
        assert (
            "opt=modredundant" in g16_fake_opt_out.route_string
            and "b3lyp empiricaldispersion=gd3bj def2svp"
            in g16_fake_opt_out.route_string
        )
        mol = Molecule.from_filepath(outputfile)
        assert np.allclose(
            mol.positions[0], [-1.04402, -2.39212, -1.17658], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
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
            "C",
            "O",
            "O",
            "H",
            "C",
            "N",
            "N",
            "N",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "C",
            "C",
            "C",
            "Cl",
            "C",
            "Cl",
            "H",
            "H",
            "O",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "C",
            "H",
            "Cl",
            "C",
            "H",
            "H",
            "H",
        ]

    def test_gaussian_scan_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2], [3, 4, 5]],
            "num_steps": 10,
            "step_size": 0.1,
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(tmpdir, "gaussian_scan_fake.com")

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "opt=modredundant" in g16_fake_opt_in.route_string
            and "b3lyp empiricaldispersion=gd3bj def2svp"
            in g16_fake_opt_in.route_string
        )
        assert (
            "B 1 2 S 10 0.1" in g16_fake_opt_in.contents
            and "A 3 4 5 S 10 0.1" in g16_fake_opt_in.contents
        )

        mol = Molecule.from_filepath(inputfile)
        assert np.allclose(
            mol.positions[0], [-1.04402, -2.39212, -1.17658], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
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
            "C",
            "O",
            "O",
            "H",
            "C",
            "N",
            "N",
            "N",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "C",
            "C",
            "C",
            "Cl",
            "C",
            "Cl",
            "H",
            "H",
            "O",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "C",
            "H",
            "Cl",
            "C",
            "H",
            "H",
            "H",
        ]

    def test_gaussian_ts_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_ts_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.ts_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianTSJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_ts",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(tmpdir, "gaussian_ts_fake.com")

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "opt=(ts,calcfc,noeigentest)" in g16_fake_opt_in.route_string
            and "b3lyp empiricaldispersion=gd3bj def2svp"
            in g16_fake_opt_in.route_string
        )
        mol = Molecule.from_filepath(inputfile)
        assert np.allclose(
            mol.positions[0], [-1.04402, -2.39212, -1.17658], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
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
            "C",
            "O",
            "O",
            "H",
            "C",
            "N",
            "N",
            "N",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "H",
            "C",
            "C",
            "C",
            "Cl",
            "C",
            "Cl",
            "H",
            "H",
            "O",
            "C",
            "H",
            "C",
            "C",
            "C",
            "H",
            "H",
            "H",
            "H",
            "C",
            "C",
            "H",
            "Cl",
            "C",
            "H",
            "H",
            "H",
        ]

    def test_gaussian_sp_input_with_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_sp_from_nhc_singlet_log_with_solvent_file,
    ):
        """Test writing simple .com input file using settings from .log file,
        including solvation."""

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)
        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_from_log_with_solvent",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_solvent_fake.com"
        )

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "scrf=(smd,solvent=toluene)" in g16_fake_opt_in.route_string
            and "b3lyp empiricaldispersion=gd3bj def2tzvp"
            in g16_fake_opt_in.route_string
        )
        mol = Molecule.from_filepath(inputfile)
        assert np.allclose(
            mol.positions[0], [0.92007, 0.547051, -0.992868], atol=1e-5
        )
        assert mol.chemical_symbols == [
            "C",
            "C",
            "C",
            "H",
            "N",
            "N",
            "N",
            "C",
            "H",
            "H",
            "C",
            "O",
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
            "I",
            "I",
            "C",
            "F",
            "F",
            "F",
        ]

    def test_gaussian_sp_with_custom_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        smd_TBME_solvent_parameters_text_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file,
    ):
        """Test writing input file from log file.
        Simply taking the Gaussian nhc_neutral_singlet.log output and write
        gaussian_sp_custom_solv.com using the settings from the .log
        file and including custom solvation parameters from file smd_TBME.
        """
        smd_TBME_tmp_path = os.path.join(tmpdir, "smd_TBME")
        copy(smd_TBME_solvent_parameters_text_file, smd_TBME_tmp_path)

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        sp_settings.solvent_model = "smd"
        sp_settings.custom_solvent = smd_TBME_tmp_path
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)

        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_from_log_with_custom_solvent",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_solvent_fake.com"
        )

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "scrf=(smd,solvent=toluene)" in g16_fake_opt_in.route_string
            and "b3lyp empiricaldispersion=gd3bj def2tzvp"
            in g16_fake_opt_in.route_string
        )
        # check the custom solvent parameters are included
        assert "stoichiometry=C5H12O" in g16_fake_opt_in.contents
        assert "epsinf=1.86704896" in g16_fake_opt_in.contents
        assert "eps=2.6" in g16_fake_opt_in.contents
        assert "HBondAcidity=0.00" in g16_fake_opt_in.contents
        assert "HBondBasicity=0.54" in g16_fake_opt_in.contents
        assert (
            "SurfaceTensionAtInterface=15.7173744" in g16_fake_opt_in.contents
        )
        assert "CarbonAromaticity=0.0" in g16_fake_opt_in.contents
        assert "ElectronegativeHalogenicity=0.00" in g16_fake_opt_in.contents

    def test_gaussian_ts_with_custom_basis_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        Ni_def2tzvp_PCHOSi_svp_text_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_file,
    ):
        custom_basis_tmp_path = os.path.join(tmpdir, "custom_basis.txt")
        copy(Ni_def2tzvp_PCHOSi_svp_text_file, custom_basis_tmp_path)

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        ts_settings.gen_genecp_file = custom_basis_tmp_path
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        ts_settings = ts_settings.merge(job_settings)

        job = GaussianTSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="gaussian_sp_from_log_with_custom_basis",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_basis_fake.com"
        )

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "b3lyp empiricaldispersion=gd3bj gen"
            in g16_fake_opt_in.route_string
        )
        # check the custom solvent parameters are included
        assert (
            "!  Def2-TZVPPD  EMSL  Basis Set Exchange Library   3/8/19 9:00 AM"
            in g16_fake_opt_in.contents
        )
        assert "Ni     0" in g16_fake_opt_in.contents
        assert (
            "! Elements                             References"
            in g16_fake_opt_in.contents
        )

    def test_gaussian_ts_with_custom_basis_using_api(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_ts_genecp_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_ts_genecp_outfile
        )
        ts_settings = ts_settings.merge(job_settings)

        # update settings for heavy elements
        ts_settings.heavy_elements = ["Pd"]
        ts_settings.heavy_elements_basis = "def2-TZVPPD"
        ts_settings.light_elements_basis = "def2-SVP"

        molecule = Molecule.from_filepath(gaussian_ts_genecp_outfile)

        job = GaussianTSJob(
            molecule=molecule,
            settings=ts_settings,
            label="gaussian_sp_from_log_with_custom_basis_from_api",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_basis_from_api_fake.com"
        )

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "b3lyp empiricaldispersion=gd3bj def2svp"
            in g16_fake_opt_in.route_string
            and "opt=(ts,calcfc,noeigentest)" in g16_fake_opt_in.route_string
        )
        # check the custom solvent parameters are included
        assert "! Basis Set Exchange" in g16_fake_opt_in.contents
        assert "!   Basis set: def2-TZVPPD" in g16_fake_opt_in.contents
        assert "Pd     0" in g16_fake_opt_in.contents
        assert "PD-ECP     3     28" in g16_fake_opt_in.contents

    def test_gaussian_modred_with_custom_basis_for_all_elements_in_structure_using_api(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        modred_genecp_inputfile,
        gaussian_jobrunner_no_scratch,
        gaussian_modred_with_custom_basis_for_all_atoms_from_api,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        modred_settings = project_settings.modred_settings()
        job_settings = GaussianJobSettings.from_filepath(
            modred_genecp_inputfile
        )
        modred_settings = modred_settings.merge(job_settings)

        modred_settings.modred = [[1, 2], [3, 4, 5]]
        modred_settings.basis = "genecp"
        modred_settings.heavy_elements = [
            "C",
            "H",
            "O",
            "N",
            "Pd",
            "P",
            "S",
            "F",
            "Br",
            "I",
        ]
        # more than all elements in the system but will be filtered to only those in the system for input preparation
        modred_settings.heavy_elements_basis = "def2-TZVPPD"
        modred_settings.light_elements_basis = None  # light element basis not specified as all use custom basis from heavy_elements_basis

        molecule = Molecule.from_filepath(modred_genecp_inputfile)
        job = GaussianModredJob(
            molecule=molecule,
            settings=modred_settings,
            label="gaussian_modred_with_custom_basis_for_all_atoms_from_api",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianModredJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(
            tmpdir,
            "gaussian_modred_with_custom_basis_for_all_atoms_from_api_fake.com",
        )

        g16_fake_opt_in = Gaussian16Input(inputfile)
        assert (
            "b3lyp empiricaldispersion=gd3bj genecp"
            in g16_fake_opt_in.route_string
            and "opt=modredundant" in g16_fake_opt_in.route_string
        )
        # check the custom solvent parameters are included
        assert "! Basis Set Exchange" in g16_fake_opt_in.contents
        assert "!   Basis set: def2-TZVPPD" in g16_fake_opt_in.contents
        assert "Pd     0" in g16_fake_opt_in.contents
        assert "PD-ECP     3     28" in g16_fake_opt_in.contents

    def test_gaussian_job_pbc_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_pbc_2d_outputfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_from_graphite_2d_pbc_log,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        opt_settings = project_settings.opt_settings()
        file_settings = GaussianJobSettings.from_filepath(
            filepath=gaussian_pbc_2d_outputfile
        )
        opt_settings = opt_settings.merge(file_settings)

        molecule = Molecule.from_filepath(gaussian_pbc_2d_outputfile)
        g16_outputfile = Gaussian16Output(gaussian_pbc_2d_outputfile)
        assert g16_outputfile.list_of_pbc_conditions == [1, 1, 0]
        assert g16_outputfile.input_translation_vectors == [
            [2.47532, 0.0, 0.0],
            [-1.21995, 2.13345, 0.0],
        ]

        job = GaussianSinglePointJob(
            molecule=molecule,
            settings=opt_settings,
            label="graphite_2d_opt_from_log",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)

        # set job.folder to the tmpdir
        job.folder = tmpdir
        job.run()  # run the job to create the output file
        assert job.is_complete()
        assert os.path.isfile(job.outputfile)
        inputfile = os.path.join(tmpdir, "graphite_2d_opt_from_log_fake.com")

        g16_fake_opt_in = Gaussian16Input(inputfile)

        assert (
            "b3lyp empiricaldispersion=gd3bj def2svp"
            in g16_fake_opt_in.route_string
            and "opt freq" in g16_fake_opt_in.route_string
        )
        assert (
            "TV       2.4755330000   -0.0060580000   -0.0000000000"
            in g16_fake_opt_in.contents
        )
