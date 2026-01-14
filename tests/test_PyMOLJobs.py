import os.path
import shutil

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol.align import PyMOLAlignJob
from chemsmart.jobs.mol.irc import PyMOLIRCMovieJob
from chemsmart.jobs.mol.mo import PyMOLMOJob
from chemsmart.jobs.mol.movie import PyMOLMovieJob
from chemsmart.jobs.mol.nci import PyMOLNCIJob
from chemsmart.jobs.mol.spin import PyMOLSpinJob
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

    def test_pymol_align_job_on_three_files(
        self,
        tmpdir,
        orca_input_nebts_reactant_xyz_file,
        gaussian_frozen_opt_inputfile,
        gaussian_singlet_opt_outfile,
        pymol_align_jobrunner,
    ):
        mol1 = Molecule.from_filepath(orca_input_nebts_reactant_xyz_file)
        mol2 = Molecule.from_filepath(gaussian_frozen_opt_inputfile)
        mol3 = Molecule.from_filepath(gaussian_singlet_opt_outfile)
        mol1.name = "R-1a_opt"
        mol2.name = "frozen_coordinates_opt"
        mol3.name = "nhc_neutral_singlet"

        job = PyMOLAlignJob(
            molecule=[mol1, mol2, mol3],
            label="R-1a_opt_and_2_molecules_align",
            jobrunner=pymol_align_jobrunner,
        )
        job.set_folder(tmpdir)
        job.run()

        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        pse_file = os.path.join(tmpdir, "R-1a_opt_and_2_molecules_align.pse")
        mol1_xyz = os.path.join(tmpdir, "R-1a_opt.xyz")
        mol2_xyz = os.path.join(tmpdir, "frozen_coordinates_opt.xyz")
        mol3_xyz = os.path.join(tmpdir, "nhc_neutral_singlet.xyz")

        assert os.path.exists(style_file)
        assert os.path.exists(pse_file)
        assert os.path.exists(mol1_xyz)
        assert os.path.exists(mol2_xyz)
        assert os.path.exists(mol3_xyz)

        out_file = os.path.join(tmpdir, "R-1a_opt_and_2_molecules_align.out")

        if os.path.exists(out_file):
            with open(out_file, "r") as f:
                content = f.read()
                assert "align frozen_coordinates_opt, R-1a_opt" in content
                assert "align nhc_neutral_singlet, R-1a_opt" in content

        for xyz_file in [mol1_xyz, mol2_xyz, mol3_xyz]:
            assert os.path.getsize(xyz_file) > 0

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

    def test_pymol_irc_movie_job_on_gaussian_log_multiple_structures(
        self,
        tmpdir,
        gaussian_singlet_opt_outfile,
        pymol_ircmovie_jobrunner,
    ):
        molecules = Molecule.from_filepath(
            gaussian_singlet_opt_outfile, index=":", return_list=True
        )

        job = PyMOLIRCMovieJob(
            molecules=molecules,
            label="nhc_neutral_singlet",
            jobrunner=pymol_ircmovie_jobrunner,
        )
        job.set_folder(tmpdir)

        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "nhc_neutral_singlet_movie.xyz")
        pse_file = os.path.join(tmpdir, "nhc_neutral_singlet_movie.pse")
        movie_file = os.path.join(tmpdir, "nhc_neutral_singlet_movie.mp4")
        assert os.path.exists(style_file)
        assert os.path.exists(xyz_file)
        assert os.path.exists(pse_file)
        assert os.path.exists(movie_file)

        molecules_check = Molecule.from_filepath(
            xyz_file, index=":", return_list=True
        )
        assert (
            len(molecules_check) == 10
        ), f"Expected 10 molecules, but got {len(molecules_check)}."

    def test_pymol_MO_job_parameters(
        self,
        tmpdir,
        gaussian_benzene_opt_outfile,
    ):

        molecules = Molecule.from_filepath(
            gaussian_benzene_opt_outfile, index="-1", return_list=True
        )

        job_homo = PyMOLMOJob(
            molecules,
            label="benzene",
            homo=True,
            lumo=False,
            number=None,
        )
        job_homo.set_folder(tmpdir)

        assert job_homo.homo is True
        assert job_homo.lumo is False
        assert job_homo.number is None
        assert job_homo.label == "benzene"
        assert job_homo.mo_basename == "benzene_HOMO"
        assert job_homo.TYPE == "pymol_mo"

        job_lumo = PyMOLMOJob(
            molecules,
            label="benzene",
            homo=False,
            lumo=True,
            number=None,
        )
        job_lumo.set_folder(tmpdir)

        assert job_lumo.homo is False
        assert job_lumo.lumo is True
        assert job_lumo.number is None
        assert job_lumo.label == "benzene"
        assert job_lumo.mo_basename == "benzene_LUMO"
        assert job_lumo.TYPE == "pymol_mo"

        job_mo5 = PyMOLMOJob(
            molecules,
            label="benzene",
            homo=False,
            lumo=False,
            number=5,
        )

        assert job_mo5.homo is False
        assert job_mo5.lumo is False
        assert job_mo5.number == 5
        assert job_mo5.label == "benzene"
        assert job_mo5.mo_basename == "benzene_MO5"

    def test_pymol_spin_job_parameters(
        self,
        tmpdir,
        gaussian_benzene_opt_outfile,
    ):
        molecules = Molecule.from_filepath(
            gaussian_benzene_opt_outfile, index="-1", return_list=True
        )

        job_spin_default = PyMOLSpinJob(
            molecules,
            label="benzene_spin",
            npts=80,
        )
        job_spin_default.set_folder(tmpdir)
        assert job_spin_default.npts == 80
        assert job_spin_default.label == "benzene_spin"
        assert job_spin_default.spin_basename == "benzene_spin_spin"
        assert job_spin_default.TYPE == "pymol_spin"
        assert not job_spin_default.is_complete()

        job_spin_custom = PyMOLSpinJob(
            molecules,
            label="benzene_spin",
            npts=100,
        )
        job_spin_custom.set_folder(tmpdir)
        assert job_spin_custom.npts == 100
        assert job_spin_custom.label == "benzene_spin"
        assert job_spin_custom.spin_basename == "benzene_spin_spin"
        assert job_spin_custom.TYPE == "pymol_spin"

        job_spin_string = PyMOLSpinJob(
            molecules,
            label="benzene_spin",
            npts="120 h",
        )
        assert job_spin_string.npts == "120 h"
        assert job_spin_string.label == "benzene_spin"
        assert job_spin_string.spin_basename == "benzene_spin_spin"
        assert job_spin_string.TYPE == "pymol_spin"

    def test_pymol_nci_job_parameters(
        self,
        tmpdir,
        gaussian_benzene_opt_outfile,
    ):
        molecules = Molecule.from_filepath(
            gaussian_benzene_opt_outfile, index="-1", return_list=True
        )

        job_nci_default = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=0.5,
            color_range=1.0,
            binary=False,
            intermediate=False,
        )
        job_nci_default.set_folder(tmpdir)
        assert job_nci_default.binary is False
        assert job_nci_default.intermediate is False
        assert job_nci_default.label == "benzene"
        assert job_nci_default.nci_basename == "benzene_nci"
        assert job_nci_default.isosurface_value == 0.5
        assert job_nci_default.color_range == 1.0
        assert job_nci_default.TYPE == "pymol_nci"
        assert not job_nci_default.is_complete()

        job_nci_binary = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=0.5,
            color_range=1.0,
            binary=True,
            intermediate=False,
        )
        job_nci_binary.set_folder(tmpdir)
        assert job_nci_binary.binary is True
        assert job_nci_binary.intermediate is False
        assert job_nci_binary.label == "benzene"
        assert job_nci_binary.nci_basename == "benzene_nci_binary"
        assert job_nci_binary.TYPE == "pymol_nci"

        job_nci_intermediate = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=0.5,
            color_range=1.0,
            binary=False,
            intermediate=True,
        )
        assert job_nci_intermediate.binary is False
        assert job_nci_intermediate.intermediate is True
        assert job_nci_intermediate.label == "benzene"
        assert job_nci_intermediate.nci_basename == "benzene_nci_intermediate"
        assert job_nci_intermediate.TYPE == "pymol_nci"

        job_nci_custom = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=0.3,
            color_range=2.0,
        )
        assert job_nci_custom.isosurface_value == 0.3
        assert job_nci_custom.color_range == 2.0
        assert job_nci_custom.label == "benzene"
        assert job_nci_custom.nci_basename == "benzene_nci"

        job_nci_combined = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=0.5,
            color_range=1.0,
            binary=True,
            intermediate=True,
        )
        assert job_nci_combined.binary is True
        assert job_nci_combined.intermediate is True
        assert "binary" in job_nci_combined.nci_basename
        assert "intermediate" in job_nci_combined.nci_basename

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
