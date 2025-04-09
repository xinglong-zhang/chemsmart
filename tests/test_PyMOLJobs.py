import os.path

from chemsmart.jobs.mol.visualize import PyMOLVisualizationJob


class TestPyMOLJobs:
    def test_pymol_visualization_job_on_gaussian_com_file(
        self,
        tmpdir,
        gaussian_opt_inputfile,
        pymol_visualization_jobrunner,
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_filename(gaussian_opt_inputfile)
        job.set_folder(tmpdir)
        job.runner = pymol_visualization_jobrunner

        # run job
        job.run(jobrunner=pymol_visualization_jobrunner)
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "model_opt_input.xyz")
        pse_file = os.path.join(tmpdir, "model_opt_input.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)
