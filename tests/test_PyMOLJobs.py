import os.path
import shutil

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol import (
    PyMOLHybridVisualizationJob,
    PyMOLMovieJob,
    PyMOLNCIJob,
    PyMOLSpinJob,
)
from chemsmart.jobs.mol.visualize import PyMOLVisualizationJob
from chemsmart.utils.cluster import is_pubchem_network_available


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

    @pytest.mark.skipif(
        not is_pubchem_network_available(),
        reason="Network to pubchem is unavailable",
    )
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

    @pytest.mark.skipif(
        not is_pubchem_network_available(),
        reason="Network to pubchem is unavailable",
    )
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

    def test_pymol_hybrid_visualization_job_on_xyz_file(
        self,
        tmpdir,
        dna_hybrid_visualized_xyz_file,
        pymol_hybrid_visualization_jobrunner,
    ):
        group1 = "503-523"
        group2 = "336, 397-412, 414-422"
        group3 = "467-495, 497-500, 502"
        group4 = "524-539"
        # set up jobs
        job = PyMOLHybridVisualizationJob.from_filename(
            dna_hybrid_visualized_xyz_file,
            jobrunner=pymol_hybrid_visualization_jobrunner,
            group1=group1,
            group2=group2,
            group3=group3,
            group4=group4,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        pse_file = os.path.join(tmpdir, "dna_hybrid.pse")
        pml_file = os.path.join(tmpdir, "hybrid_visualization.pml")
        group_selection_commands = [
            "pymol_style all\n",
            "unset stick_color, all\n",
            "hide everything, all\n",
            "show sticks, all\n",
            "set_color light_C, [0.8, 0.8, 0.9]\n",
            "set_color light_N, [0.6, 0.8, 1.0]\n",
            "set_color light_O, [1.0, 0.7, 0.7]\n",
            "set_color light_P, [1.0, 0.85, 0.6]\n",
            "color light_C, elem C\n",
            "color light_P, elem P\n",
            "color light_O, elem O\n",
            "color light_N, elem N\n",
            "select group1,  id 503-523\n",
            "util.cbap group1\n",
            "select group2,  id 336 or id 397-412 or id 414-422\n",
            "util.cbac group2\n",
            "select group3,  id 467-495 or id 497-500 or id 502\n",
            "util.cbay group3\n",
            "select group4,  id 524-539\n",
            "util.cbag group4\n",
            "set stick_radius, 0.25, (group1 or group2 or group3 or group4)\n",
        ]
        with open(pml_file, "r") as f:
            content = f.readlines()
            print(content)
            for i in group_selection_commands:
                assert i in content
        assert os.path.exists(style_file)
        assert os.path.exists(pse_file)

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
        assert job.job_basename == "model_opt_input_movie"
        assert job.outputfile == os.path.join(
            tmpdir, "model_opt_input_movie.pse"
        )

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "model_opt_input_movie.xyz")
        pse_file = os.path.join(tmpdir, "model_opt_input_movie.pse")
        movie_file = os.path.join(tmpdir, "model_opt_input_movie.mp4")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)
        assert os.path.exists(movie_file)

    def test_pymol_nci_job_on_gaussian_com_file(
        self,
        tmpdir,
        gaussian_opt_inputfile,
        pymol_movie_jobrunner,
    ):
        # set up jobs
        job = PyMOLNCIJob.from_filename(
            gaussian_opt_inputfile,
            jobrunner=pymol_movie_jobrunner,
            isosurface_value=0.5,
            color_range=1.2,
        )
        job.set_folder(tmpdir)
        assert job.job_basename == "model_opt_input_nci"
        assert job.outputfile == os.path.join(
            tmpdir, "model_opt_input_nci.pse"
        )

    def test_pymol_spin_job_on_gaussian_com_file(
        self,
        tmpdir,
        gaussian_opt_inputfile,
        pymol_movie_jobrunner,
    ):
        # set up jobs
        job = PyMOLSpinJob.from_filename(
            gaussian_opt_inputfile,
            jobrunner=pymol_movie_jobrunner,
            isosurface_value=0.5,
            color_range=1.2,
        )
        job.set_folder(tmpdir)
        assert job.job_basename == "model_opt_input_spin"
        assert job.outputfile == os.path.join(
            tmpdir, "model_opt_input_spin.pse"
        )
