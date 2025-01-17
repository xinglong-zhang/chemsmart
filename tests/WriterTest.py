import os
from filecmp import cmp
from shutil import copy

import pytest
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.settings import read_molecular_job_yaml
from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.server import Server
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


class TestGaussianInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        pbs_server,
        jobrunner_no_scratch,
    ):
        # get project settings
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
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        assert os.path.isfile(os.path.join(tmpdir, "gaussian_opt.com"))
        print(os.path.join(tmpdir, "gaussian_opt.com"))
