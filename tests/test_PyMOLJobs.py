# class TestPyMOLJobs:
#     def test_pymol_visualization_job_on_gaussian_com_file(
#         self,
#         tmpdir,
#         gaussian_opt_inputfile,
#         gaussian_yaml_settings_gas_solv_project_name,
#         jobrunner_no_scratch,
#         gaussian_written_opt_file,
#     ):
#
#         job = PyMOLVisualizationJob.from_filename(gaussian_opt_inputfile)
#
#         # get project settings
#         project_settings = GaussianProjectSettings.from_project(
#             gaussian_yaml_settings_gas_solv_project_name
#         )
#         settings = project_settings.opt_settings()
#         settings.charge = 0
#         settings.multiplicity = 1
#         job = GaussianOptJob.from_filename(
#             filename=single_molecule_xyz_file,
#             settings=settings,
#             label="gaussian_opt",
#         )
#         assert isinstance(job, GaussianOptJob)
#         g16_writer = GaussianInputWriter(
#             job=job, jobrunner=jobrunner_no_scratch
#         )
#
#         # write input file
#         g16_writer.write(target_directory=tmpdir)
#         g16_file = os.path.join(tmpdir, "gaussian_opt.com")
#         assert os.path.isfile(g16_file)
#         assert cmp(
#             g16_file, gaussian_written_opt_file, shallow=False
#         )  # writes input file as expected
#
#         # job run will result in the job being run and the output file copied back to run folder
#         # job.run(jobrunner=jobrunner_no_scratch)
#         # assert job.is_complete()
