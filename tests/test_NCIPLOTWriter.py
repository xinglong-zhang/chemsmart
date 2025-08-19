# import os
# from filecmp import cmp
# from shutil import copy
#
# from chemsmart.io.gaussian.output import Gaussian16Output
# from chemsmart.io.molecules.structure import Molecule
# from chemsmart.jobs.gaussian import (
#     GaussianModredJob,
#     GaussianOptJob,
#     GaussianScanJob,
#     GaussianSinglePointJob,
#     GaussianTSJob,
# )
# from chemsmart.jobs.gaussian.settings import GaussianJobSettings
# from chemsmart.jobs.gaussian.writer import GaussianInputWriter
# from chemsmart.jobs.nciplot import NCIPLOTJob
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings

# from chemsmart.settings.gaussian import GaussianProjectSettings
# from chemsmart.utils.utils import cmp_with_ignore


class TestNCIPLOTInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file,
    ):
        job_settings = NCIPLOTJobSettings(
            filenames=single_molecule_xyz_file,
            label="gaussian_opt",
            rthres=10.0,
            ligand_file_number=1,
            ligand_radius=1.0,
            radius_positions=[1, 2, 3],
            radius_r=1.5,
            intercut1=0.5,
            intercut2=0.5,
            increments="0.1,0.1,0.1",
            fragments=1,
            cutoff_density_dat=0.01,
            cutoff_rdg_dat=0.01,
            cutoff_density_cube=0.01,
            cutoff_rdg_cube=0.01,
            dgrid=0.1,
            integrate=True,
            ranges=[0.0, 1.0, 2.0],
        )
        print(job_settings)
        #
        # # get project settings
        # project_settings = GaussianProjectSettings.from_project(
        #     gaussian_yaml_settings_gas_solv_project_name
        # )
        # settings = project_settings.opt_settings()
        # settings.charge = 0
        # settings.multiplicity = 1
        # job = NCIPLOTJob(
        #     filenames=single_molecule_xyz_file,
        #     settings=settings,
        #     label="gaussian_opt",
        #     jobrunner=gaussian_jobrunner_no_scratch,
        # )
        #
        # assert isinstance(job, GaussianOptJob)
        # g16_writer = GaussianInputWriter(job=job)
        #
        # # write input file
        # g16_writer.write(target_directory=tmpdir)
        # g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        # assert os.path.isfile(g16_file)
        # assert cmp(
        #     g16_file, gaussian_written_opt_file, shallow=False
        # )  # writes input file as expected
        #
        # # job run will result in the job being run and the output file copied back to run folder
        # # job.run()
        # # assert job.is_complete()
