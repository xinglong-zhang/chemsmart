import os.path
import shutil

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol import PyMOLHybridVisualizationJob
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
            label="phenyldioxazolone",
            # pymol label should avoid "," which is a separator for commands,
            # such as e.g., 3-Phenyl-1,4,2-dioxazol-5-one
            jobrunner=pymol_visualization_jobrunner,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        xyz_file = os.path.join(tmpdir, "phenyldioxazolone.xyz")
        pse_file = os.path.join(tmpdir, "phenyldioxazolone.pse")
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
        group5 = "541-550"
        groups = [group1, group2, group3, group4, group5]
        # set up jobs
        job = PyMOLHybridVisualizationJob.from_filename(
            dna_hybrid_visualized_xyz_file,
            jobrunner=pymol_hybrid_visualization_jobrunner,
            groups=groups,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        pse_file = os.path.join(tmpdir, "dna_hybrid_hybrid_visualization.pse")
        pml_file = os.path.join(
            tmpdir, f"{os.path.basename(tmpdir)}_hybrid_visualization.pml"
        )
        group_selection_commands = [
            "pymol_style all\n",
            "unset stick_color, all\n",
            "hide everything, all\n",
            "show sticks, all\n",
            "set_color light_C, [0.8, 0.8, 0.9]\n",
            "set_color light_N, [0.6, 0.8, 1.0]\n",
            "set_color light_O, [1.0, 0.7, 0.7]\n",
            "set_color light_P, [1.0, 0.85, 0.6]\n",
            "set_color light_S, [1.0, 0.7, 0.7]\n",
            "color light_C, elem C\n",
            "color light_P, elem P\n",
            "color light_O, elem O\n",
            "color light_N, elem N\n",
            "color light_S, elem S\n",
            "select group1,  id 503-523\n",
            "select group2,  id 336 or id 397-412 or id 414-422\n",
            "select group3,  id 467-495 or id 497-500 or id 502\n",
            "select group4,  id 524-539\n",
            "select group5,  id 541-550\n",
            "util.cbap group1\n",
            "util.cbac group2\n",
            "util.cbay group3\n",
            "util.cbag group4\n",
            "util.cbam group5\n",
            "set stick_transparency, 0, all\n",
            "set stick_radius, 0.25, (group1 or group2 or group3 or group4 or group5)\n",
            "show surface, all\n",
            "set surface_color, grey, all\n",
            "set transparency, 0.7, all\n",
        ]
        with open(pml_file, "r") as f:
            content = f.readlines()
            for i in group_selection_commands:
                assert i in content
        assert os.path.exists(style_file)
        assert os.path.exists(pse_file)

    def test_pymol_hybrid_visualization_job_with_redundant_colors_on_xyz_file(
        self,
        tmpdir,
        dna_hybrid_visualized_xyz_file,
        pymol_hybrid_visualization_jobrunner,
    ):
        group1 = "503-523"
        group2 = "336, 397-412, 414-422"
        group3 = "467-495, 497-500, 502"
        color1 = "cbap"
        color2 = "cbak"
        color3 = "cbam"
        color4 = "cbay"
        groups = [group1, group2, group3]
        colors = [color1, color2, color3, color4]
        # set up jobs
        job = PyMOLHybridVisualizationJob.from_filename(
            dna_hybrid_visualized_xyz_file,
            jobrunner=pymol_hybrid_visualization_jobrunner,
            groups=groups,
            colors=colors,
        )
        job.set_folder(tmpdir)

        # run job
        job.run()
        assert job.is_complete()
        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        pse_file = os.path.join(tmpdir, "dna_hybrid_hybrid_visualization.pse")
        pml_file = os.path.join(
            tmpdir, f"{os.path.basename(tmpdir)}_hybrid_visualization.pml"
        )
        group_selection_commands = [
            "select group1,  id 503-523\n",
            "select group2,  id 336 or id 397-412 or id 414-422\n",
            "select group3,  id 467-495 or id 497-500 or id 502\n",
            "util.cbap group1\n",
            "util.cbak group2\n",
            "util.cbam group3\n",
            "set stick_transparency, 0, all\n",
            "set stick_radius, 0.25, (group1 or group2 or group3)\n",
            "show surface, all\n",
            "set surface_color, grey, all\n",
            "set transparency, 0.7, all\n",
        ]
        with open(pml_file, "r") as f:
            content = f.readlines()
            for i in group_selection_commands:
                assert i in content
        assert os.path.exists(style_file)
        assert os.path.exists(pse_file)

    def test_pymol_hybrid_visualization_job_custom_light_colors_on_xyz_file(
        self,
        tmpdir,
        dna_hybrid_visualized_xyz_file,
        pymol_hybrid_visualization_jobrunner,
    ):
        # verify that custom light colors provided to the job are written to the pml
        group1 = "503-523"
        group2 = "336, 397-412, 414-422"
        group3 = "467-495, 497-500, 502"
        groups = [group1, group2, group3]

        # custom RGB values for light colors
        light_colors = {
            "C": [0.1, 0.2, 0.3],
            "N": [0.2, 0.3, 0.4],
            "O": [0.3, 0.4, 0.5],
            "P": [0.6, 0.7, 0.8],
            "S": [0.9, 0.8, 0.7],
        }
        new_color_carbon = light_colors["C"]
        new_color_nitrogen = light_colors["N"]
        new_color_oxygen = light_colors["O"]
        new_color_phosphorus = light_colors["P"]
        new_color_sulfur = light_colors["S"]

        job = PyMOLHybridVisualizationJob.from_filename(
            dna_hybrid_visualized_xyz_file,
            jobrunner=pymol_hybrid_visualization_jobrunner,
            groups=groups,
            new_color_carbon=new_color_carbon,
            new_color_nitrogen=new_color_nitrogen,
            new_color_oxygen=new_color_oxygen,
            new_color_phosphorus=new_color_phosphorus,
            new_color_sulfur=new_color_sulfur,
        )
        job.set_folder(tmpdir)

        job.run()
        assert job.is_complete()

        pml_file = os.path.join(
            tmpdir, f"{os.path.basename(tmpdir)}_hybrid_visualization.pml"
        )
        group_selection_commands = [
            "set_color light_C, [0.1, 0.2, 0.3]\n",
            "set_color light_N, [0.2, 0.3, 0.4]\n",
            "set_color light_O, [0.3, 0.4, 0.5]\n",
            "set_color light_P, [0.6, 0.7, 0.8]\n",
            "set_color light_S, [0.9, 0.8, 0.7]\n",
        ]

        with open(pml_file, "r") as f:
            content = f.readlines()
            for line in group_selection_commands:
                assert line in content

        style_file = os.path.join(tmpdir, "zhang_group_pymol_style.py")
        pse_file = os.path.join(tmpdir, "dna_hybrid_hybrid_visualization.pse")
        assert os.path.exists(style_file)
        assert os.path.exists(pse_file)

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
