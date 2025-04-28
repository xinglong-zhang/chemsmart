import os.path
import shutil

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol import PyMOLMovieJob
from chemsmart.jobs.mol.visualize import PyMOLVisualizationJob


@pytest.fixture(scope="session")
def skip_if_no_pymol():
    if shutil.which("pymol") is None:
        pytest.skip("PyMOL not installed")


@pytest.mark.usefixtures("skip_if_no_pymol")
class TestPyMOLJobs:
    def test_pymol_visualization_job_on_gaussian_com_file(
        self,
        tmpdir,
        gaussian_opt_inputfile,
        pymol_visualization_jobrunner,
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_filename(
            gaussian_opt_inputfile, jobrunner=pymol_visualization_jobrunner
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "model_opt_input.xyz")
        pse_file = os.path.join(tmpdir, "model_opt_input.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

    def test_pymol_visualization_job_on_gaussian_log_file(
        self,
        gaussian_singlet_opt_outfile,
        tmpdir,
        pymol_visualization_jobrunner,
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_filename(
            gaussian_singlet_opt_outfile,
            jobrunner=pymol_visualization_jobrunner,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "nhc_neutral_singlet.xyz")
        pse_file = os.path.join(tmpdir, "nhc_neutral_singlet.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

        molecules = Molecule.from_filepath(
            xyz_file, index=":", return_list=True
        )
        assert (
            len(molecules) == 1
        ), f"Expected 1 molecule, but got {len(molecules)}."

    def test_pymol_visualization_job_on_orca_inp_file(
        self, tmpdir, water_opt_input_path, pymol_visualization_jobrunner
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_filename(
            water_opt_input_path, jobrunner=pymol_visualization_jobrunner
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "water_opt.xyz")
        pse_file = os.path.join(tmpdir, "water_opt.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

    def test_pymol_visualization_job_on_orca_out_file(
        self, tmpdir, water_output_gas_path, pymol_visualization_jobrunner
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_filename(
            water_output_gas_path, jobrunner=pymol_visualization_jobrunner
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "water_opt.xyz")
        pse_file = os.path.join(tmpdir, "water_opt.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

    def test_pymol_visualization_job_on_pubchem_id(
        self, tmpdir, pymol_visualization_jobrunner
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_pubchem(
            "8028", label="thf", jobrunner=pymol_visualization_jobrunner
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "thf.xyz")
        pse_file = os.path.join(tmpdir, "thf.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

    def test_pymol_visualization_job_on_smiles(
        self, tmpdir, pymol_visualization_jobrunner
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_pubchem(
            "C1=CC=C(C=C1)C2=NOC(=O)O2",
            label="3-Phenyl-1,4,2-dioxazol-5-one",
            jobrunner=pymol_visualization_jobrunner,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "3-Phenyl-1,4,2-dioxazol-5-one.xyz")
        pse_file = os.path.join(tmpdir, "3-Phenyl-1,4,2-dioxazol-5-one.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

    def test_pymol_visualization_job_on_gaussian_log_multiple_structures(
        self,
        tmpdir,
        gaussian_singlet_opt_outfile,
        pymol_visualization_jobrunner,
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_filename(
            gaussian_singlet_opt_outfile,
            index=":",
            jobrunner=pymol_visualization_jobrunner,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "nhc_neutral_singlet.xyz")
        pse_file = os.path.join(tmpdir, "nhc_neutral_singlet.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)

        molecules = Molecule.from_filepath(
            xyz_file, index=":", return_list=True
        )
        assert (
            len(molecules) == 10
        ), f"Expected 1 molecule, but got {len(molecules)}."

    def test_pymol_movie_job_on_gaussian_com_file(
        self,
        tmpdir,
        gaussian_opt_inputfile,
        pymol_movie_jobrunner,
    ):
        # set up jobs
        job = PyMOLMovieJob.from_filename(
            gaussian_opt_inputfile, jobrunner=pymol_movie_jobrunner
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "model_opt_input.xyz")
        pse_file = os.path.join(tmpdir, "model_opt_input.pse")
        movie_file = os.path.join(tmpdir, "model_opt_input.mp4")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)
        assert os.path.exists(movie_file)
