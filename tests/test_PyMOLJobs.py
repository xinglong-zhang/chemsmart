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
        print(tmpdir)
