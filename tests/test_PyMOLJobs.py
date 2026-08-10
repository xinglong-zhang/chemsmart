import io
import os.path
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol import PyMOLHybridVisualizationJob
from chemsmart.jobs.mol.align import PyMOLAlignJob
from chemsmart.jobs.mol.irc import PyMOLIRCMovieJob
from chemsmart.jobs.mol.mo import PyMOLMOJob
from chemsmart.jobs.mol.movie import PyMOLMovieJob
from chemsmart.jobs.mol.nci import PyMOLNCIJob
from chemsmart.jobs.mol.runner import (
    PYMOL_SCIENTIFIC_STYLE_COMMANDS,
    PYMOL_VISUALIZE_STYLE_CLI_CHOICES,
    PyMOLAlignJobRunner,
    PyMOLHybridVisualizationJobRunner,
    PyMOLJobRunner,
    PyMOLNCIJobRunner,
    PyMOLScientificStyleVisualizationJobRunner,
    PyMOLSpinJobRunner,
    is_pymol_derived_style,
    normalize_pymol_style,
)
from chemsmart.jobs.mol.spin import PyMOLSpinJob
from chemsmart.jobs.mol.templates import zhang_group_scientific_styles
from chemsmart.jobs.mol.templates.zhang_group_scientific_styles import (
    SCIENTIFIC_STYLE_CLASSES,
    ComicMetallicStyle,
    GlossyStyle,
    NeonCoordinationCoreStyle,
    QuasiChemDrawBoldStyle,
    ScientificStyle,
    SoftCartoonStyle,
    StericSurfaceStyle,
    metal_pymol_selection,
    pymol_elem_selection,
)
from chemsmart.jobs.mol.visualize import (
    PyMOLScientificStyleVisualizationJob,
    PyMOLVisualizationJob,
)
from chemsmart.utils.cluster import (
    is_pubchem_api_available,
    is_pubchem_network_available,
)
from chemsmart.utils.utils import quote_path


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
        not is_pubchem_network_available() or not is_pubchem_api_available(),
        reason="Network/API to pubchem is unavailable",
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
        not is_pubchem_network_available() or not is_pubchem_api_available(),
        reason="Network/API to pubchem is unavailable",
    )
    def test_pymol_visualization_job_on_smiles(
        self, tmpdir, pymol_visualization_jobrunner
    ):
        # set up jobs
        job = PyMOLVisualizationJob.from_pubchem(
            "C1=CC=C(C=C1)C2=NOC(=O)O2",
            label="phenyldioxazolone",
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
        pml_file = os.path.join(tmpdir, "dna_hybrid_hybrid_visualization.pml")
        group_selection_commands = [
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
            "select group1, id 503-523\n",
            "select group2, id 336 or id 397-412 or id 414-422\n",
            "select group3, id 467-495 or id 497-500 or id 502\n",
            "select group4, id 524-539\n",
            "select group5, id 541-550\n",
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
        pml_file = os.path.join(tmpdir, "dna_hybrid_hybrid_visualization.pml")
        group_selection_commands = [
            "select group1, id 503-523\n",
            "select group2, id 336 or id 397-412 or id 414-422\n",
            "select group3, id 467-495 or id 497-500 or id 502\n",
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
        # verify that custom light colors provided
        # to the job are written to the pml
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

        pml_file = os.path.join(tmpdir, "dna_hybrid_hybrid_visualization.pml")
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

        # None isosurface_value/color_range fall back to their defaults
        job_nci_none_defaults = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=None,
            color_range=None,
        )
        assert job_nci_none_defaults.isosurface_value == 0.5
        assert job_nci_none_defaults.color_range == 1.0

        # explicit nci_basename is used as-is (not derived from label)
        job_nci_explicit_basename = PyMOLNCIJob(
            molecules,
            label="benzene",
            isosurface_value=0.5,
            color_range=1.0,
            nci_basename="custom_name",
        )
        assert job_nci_explicit_basename.nci_basename == "custom_name"
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


class TestPyMOLCLIFolderOptions:
    """Folder options (``-d``/``-t`` and ``-d``/``-p``) in the ``mol`` CLI."""

    def test_directory_filetype_options_accepted(
        self, tmp_path, invoke_mol_with_visualize
    ):
        """``mol -d dir -t log visualize`` is accepted and populates ``ctx.obj``."""
        ctx_obj = {}
        result = invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-t", "log"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        assert ctx_obj.get("directory") == str(tmp_path)
        assert ctx_obj.get("filetype") == "log"

    def test_directory_program_options_accepted(
        self, tmp_path, invoke_mol_with_visualize
    ):
        """``mol -d dir -p gaussian visualize`` is accepted and populates ``ctx.obj``."""
        ctx_obj = {}
        result = invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-p", "gaussian"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        assert ctx_obj.get("directory") == str(tmp_path)
        assert ctx_obj.get("program") == "gaussian"

    def test_directory_filetype_label_auto_generated(
        self, tmp_path, invoke_mol_with_visualize
    ):
        """When ``-d``/``-t`` used, auto-generated label includes dir name."""
        ctx_obj = {}
        result = invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-t", "log"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        label = ctx_obj.get("label", "")
        dir_name = os.path.basename(os.path.abspath(str(tmp_path)))
        assert dir_name in label
        assert "log" in label

    def test_directory_program_label_auto_generated(
        self, tmp_path, invoke_mol_with_visualize
    ):
        """When ``-d``/``-p`` used, auto-generated label includes program name."""
        ctx_obj = {}
        result = invoke_mol_with_visualize(
            ["-d", str(tmp_path), "-p", "gaussian"],
            ctx_obj=ctx_obj,
        )
        assert "No such option" not in result.output, result.output
        label = ctx_obj.get("label", "")
        assert "gaussian" in label


class TestPyMOLFileProcessingUsesSourceFilename:
    def test_spin_cli_custom_label_uses_source_basename_and_exact_output_name(
        self, gaussian_benzene_opt_outfile, invoke_mol_cli
    ):
        from unittest.mock import patch

        custom_label = "new_name_new_spin_isovalue"
        with patch("chemsmart.jobs.mol.spin.PyMOLSpinJob") as mock_spin_job:
            result = invoke_mol_cli(
                [
                    "-f",
                    gaussian_benzene_opt_outfile,
                    "-l",
                    custom_label,
                    "spin",
                    "-i",
                    "0.1",
                ]
            )

        assert result.exit_code == 0, result.output
        _, kwargs = mock_spin_job.call_args
        assert kwargs["source_basename"] == "benzene"
        assert kwargs["label"] == custom_label
        assert kwargs["spin_basename"] == custom_label

    def test_generate_fchk_uses_source_basename_not_label(
        self,
        tmpdir,
        gaussian_benzene_opt_outfile,
        pymol_mo_jobrunner,
        monkeypatch,
    ):
        molecules = Molecule.from_filepath(
            gaussian_benzene_opt_outfile, index="-1", return_list=True
        )
        job = PyMOLMOJob(
            molecules,
            label="custom_label",
            source_basename="benzene_opt",
            homo=True,
        )
        job.set_folder(tmpdir)

        with open(os.path.join(tmpdir, "benzene_opt.chk"), "w"):
            pass

        commands = []
        monkeypatch.setattr(
            pymol_mo_jobrunner,
            "_get_gaussian_executable",
            lambda _job: "/gaussian",
        )
        monkeypatch.setattr(
            "chemsmart.jobs.mol.runner.run_command",
            lambda cmd: commands.append(cmd),
        )

        pymol_mo_jobrunner._generate_fchk_file(job)

        assert commands == ["/gaussian/formchk benzene_opt.chk"]

    def test_spin_cubegen_uses_source_basename_fchk(
        self, tmpdir, gaussian_benzene_opt_outfile, pbs_server, monkeypatch
    ):
        molecules = Molecule.from_filepath(
            gaussian_benzene_opt_outfile, index="-1", return_list=True
        )
        job = PyMOLSpinJob(
            molecules,
            label="spin_label",
            source_basename="benzene_opt",
        )
        job.set_folder(tmpdir)
        runner = PyMOLSpinJobRunner(server=pbs_server, scratch=False)

        commands = []
        monkeypatch.setattr(
            runner,
            "_get_gaussian_executable",
            lambda _job: "/gaussian",
        )
        monkeypatch.setattr(
            "chemsmart.jobs.mol.runner.run_command",
            lambda cmd: commands.append(cmd),
        )

        runner._generate_spin_cube_file(job)

        assert commands == [
            f"/gaussian/cubegen 0 spin benzene_opt.fchk spin_label_spin.cube {job.npts}"
        ]

    def test_nci_uses_source_basename_for_cube_loading_and_command(
        self, tmpdir, gaussian_benzene_opt_outfile, pbs_server
    ):
        molecules = Molecule.from_filepath(
            gaussian_benzene_opt_outfile, index="-1", return_list=True
        )
        job = PyMOLNCIJob(
            molecules,
            label="renamed_label",
            source_basename="benzene_opt",
            isosurface_value=0.5,
            color_range=1.0,
        )
        job.set_folder(tmpdir)
        runner = PyMOLNCIJobRunner(server=pbs_server, scratch=False)

        dens_file = os.path.join(tmpdir, "benzene_opt-dens.cube")
        grad_file = os.path.join(tmpdir, "benzene_opt-grad.cube")
        with open(dens_file, "w"):
            pass
        with open(grad_file, "w"):
            pass

        command = runner._load_cube_files(job, "cmd")
        command = runner._run_nci_command(job, command)

        assert f"load {quote_path(dens_file)}" in command
        assert f"load {quote_path(grad_file)}" in command
        assert "; nci benzene_opt" in command


class TestPyMOLStyleCommands:
    label_1_mer = "1-mer"
    coordination_bonds_1_mer = [
        [1, 2],
        [1, 5],
        [1, 36],
        [1, 3],
        [1, 15],
        [1, 8],
    ]

    def test_format_pymol_style_command_is_independent_of_coordinates(self):
        """``-c`` is handled via distance/angle labels, not style args."""
        for style, expected in (
            (
                "comic",
                f"comic {self.label_1_mer}",
            ),
            ("soft_cartoon", f"soft_cartoon {self.label_1_mer}"),
            (
                "neon_coordination_core",
                f"neon_coordination_core {self.label_1_mer}",
            ),
            (
                "editorial_minimal",
                f"editorial_minimal {self.label_1_mer}",
            ),
            ("soft_ceramic", f"soft_ceramic {self.label_1_mer}"),
            ("matte_clay", f"matte_clay {self.label_1_mer}"),
            (
                "glossy",
                f"glossy {self.label_1_mer}",
            ),
        ):
            for coordinates in (None, self.coordination_bonds_1_mer):
                job = SimpleNamespace(style=style, coordinates=coordinates)
                command = PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                    job, self.label_1_mer
                )
                assert command == expected

    def test_setup_style_cylview_flat(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(style="cylview-flat", label=self.label_1_mer)
        command = runner._setup_style(job, "pymol cmd")
        assert command.endswith(f' -d "cylview_flat_style {self.label_1_mer}')

        job_normalized = SimpleNamespace(
            style="cylview_flat", label=self.label_1_mer
        )
        command_normalized = runner._setup_style(job_normalized, "pymol cmd")
        assert command_normalized.endswith(
            f' -d "cylview_flat_style {self.label_1_mer}'
        )

        assert normalize_pymol_style("cylview-flat") == "cylview_flat"

    def test_align_setup_style_cylview_flat(self):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        job = SimpleNamespace(
            style="cylview-flat",
            mol_names=[self.label_1_mer, "2-mer"],
        )
        command = runner._setup_style(job, "pymol cmd")
        assert (
            ' -d "cylview_flat_style 1-mer; cylview_flat_style 2-mer'
            in command
        )

    def test_cylview_flat_visualization_job_label_suffix(self):
        """cylview-flat PSE/xyz names must not collide with cylview."""
        job_flat = PyMOLVisualizationJob(
            molecule=None,
            label="mol",
            style="cylview-flat",
        )
        assert job_flat.style == "cylview_flat"
        assert job_flat.label == "mol_cylview_flat_visualization"
        assert job_flat.job_basename == "mol_cylview_flat_visualization"

        job_cylview = PyMOLVisualizationJob(
            molecule=None,
            label="mol",
            style="cylview",
        )
        assert job_cylview.label == "mol"
        assert job_cylview.job_basename == "mol"

        job_pymol = PyMOLVisualizationJob(
            molecule=None,
            label="mol",
            style="pymol",
        )
        assert job_pymol.label == "mol"


class TestPyMOLJobRunnerBaseHelpers:
    """Direct tests for PyMOLJobRunner's private helper methods, using
    a bare instance (__new__ bypass) plus SimpleNamespace fake jobs
    since these do their own file I/O/regex work independent of a
    fully constructed job/server pipeline."""

    def test_is_pymol_derived_style_none_is_false(self):
        assert is_pymol_derived_style(None) is False

    def test_scratch_defaults_to_class_scratch_when_none(self, pbs_server):
        runner = PyMOLJobRunner(server=pbs_server)
        assert runner.scratch is PyMOLJobRunner.SCRATCH

    def test_executable_raises_when_pymol_not_on_path(self, mocker):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        mocker.patch(
            "chemsmart.jobs.mol.runner.shutil.which", return_value=None
        )
        with pytest.raises(FileNotFoundError, match="not found in PATH"):
            runner.executable

    def test_generate_visualization_style_script_overwrites_existing(
        self, tmp_path
    ):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(folder=str(tmp_path), label="testjob")
        dest = tmp_path / "zhang_group_pymol_style.py"
        dest.write_text("# existing stub\n")

        result = runner._generate_visualization_style_script(job)

        assert result == str(dest)
        assert os.path.exists(result)
        # isosurface_value/color_range default to None when absent
        assert job.isosurface_value is None
        assert job.color_range is None
        # content was overwritten (no longer the stub)
        assert dest.read_text() != "# existing stub\n"

    def test_generate_visualization_style_script_modifies_when_isosurface_set(
        self, tmp_path
    ):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(
            folder=str(tmp_path),
            label="testjob",
            isosurface_value=0.8,
            color_range=None,
        )
        result = runner._generate_visualization_style_script(job)
        assert "isosurface=0.8" in open(result).read()

    def test_modify_job_pymol_script_raises_when_missing(self, tmp_path):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(
            folder=str(tmp_path),
            label="testjob",
            isosurface_value=None,
            color_range=None,
        )
        with pytest.raises(FileNotFoundError, match="does not exist"):
            runner._modify_job_pymol_script(job)

    def test_modify_job_pymol_script_updates_isosurface_and_color_range(
        self, tmp_path
    ):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        src = runner.pymol_templates_path / "zhang_group_pymol_style.py"
        dest = tmp_path / "zhang_group_pymol_style.py"
        shutil.copy(src, dest)
        job = SimpleNamespace(
            folder=str(tmp_path),
            label="testjob",
            isosurface_value=0.8,
            color_range=2.0,
        )
        result = runner._modify_job_pymol_script(job, str(dest))
        content = open(result).read()
        assert "isosurface=0.8" in content
        assert "range=2.0" in content

    def test_modify_job_pymol_script_no_change_skips_rewrite(self, tmp_path):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        src = runner.pymol_templates_path / "zhang_group_pymol_style.py"
        dest = tmp_path / "zhang_group_pymol_style.py"
        shutil.copy(src, dest)
        original_mtime = os.path.getmtime(dest)
        job = SimpleNamespace(
            folder=str(tmp_path),
            label="testjob",
            isosurface_value=None,
            color_range=None,
        )
        result = runner._modify_job_pymol_script(job, str(dest))
        assert result == str(dest)
        assert os.path.getmtime(dest) == original_mtime

    def test_get_gaussian_executable(self, mocker):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        runner.server = SimpleNamespace(name="local")
        mock_exe_cls = mocker.patch(
            "chemsmart.jobs.mol.runner.GaussianExecutable"
        )
        mock_exe = mocker.MagicMock()
        mock_exe.executable_folder = "/opt/g16"
        mock_exe_cls.from_servername.return_value = mock_exe

        result = runner._get_gaussian_executable(SimpleNamespace())

        assert result == "/opt/g16"
        mock_exe_cls.from_servername.assert_called_once_with("local")

    def test_generate_fchk_file_raises_when_neither_file_exists(
        self, tmp_path
    ):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(folder=str(tmp_path), source_basename="mol")
        with pytest.raises(FileNotFoundError, match="is required"):
            runner._generate_fchk_file(job)

    def test_generate_fchk_file_skips_when_fchk_already_exists(
        self, tmp_path, mocker
    ):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        (tmp_path / "mol.fchk").write_text("stub")
        job = SimpleNamespace(folder=str(tmp_path), source_basename="mol")
        mocker.patch.object(
            PyMOLJobRunner,
            "_get_gaussian_executable",
            return_value="/opt/g16",
        )
        mock_run = mocker.patch("chemsmart.jobs.mol.runner.run_command")

        runner._generate_fchk_file(job)

        mock_run.assert_not_called()

    def test_write_input_list_with_non_molecule_raises(self, tmp_path):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(
            inputfile=str(tmp_path / "mol.xyz"),
            molecule=["not-a-molecule"],
        )
        with pytest.raises(ValueError, match="not of Molecule type"):
            runner._write_input(job)

    def test_write_input_non_list_non_molecule_raises(self, tmp_path):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(
            inputfile=str(tmp_path / "mol2.xyz"),
            molecule="not-a-molecule-or-list",
        )
        with pytest.raises(ValueError, match="not of Molecule type"):
            runner._write_input(job)

    def test_get_visualization_command_skips_r_flag_when_style_missing(
        self, mocker
    ):
        mocker.patch.object(
            PyMOLJobRunner,
            "executable",
            new_callable=mocker.PropertyMock,
            return_value="/usr/bin/pymol",
        )
        mocker.patch.object(
            PyMOLJobRunner,
            "_generate_visualization_style_script",
            return_value="/nonexistent/style.py",
        )
        job = SimpleNamespace(
            inputfile="/tmp/mol.xyz",
            pymol_script=None,
            label="mol",
            quiet_mode=False,
            command_line_only=False,
        )
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        command = runner._get_visualization_command(job)
        assert command == "/usr/bin/pymol /tmp/mol.xyz"

    def test_get_visualization_command_uses_existing_user_script(self, mocker):
        mocker.patch.object(
            PyMOLJobRunner,
            "executable",
            new_callable=mocker.PropertyMock,
            return_value="/usr/bin/pymol",
        )
        job = SimpleNamespace(
            inputfile="/tmp/mol.xyz",
            pymol_script="chemsmart/jobs/mol/runner.py",
            label="mol",
            quiet_mode=True,
            command_line_only=True,
        )
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        command = runner._get_visualization_command(job)
        assert " -r chemsmart/jobs/mol/runner.py" in command
        assert " -q" in command
        assert " -c" in command

    def test_get_visualization_command_missing_user_script_asserts(
        self, mocker
    ):
        mocker.patch.object(
            PyMOLJobRunner,
            "executable",
            new_callable=mocker.PropertyMock,
            return_value="/usr/bin/pymol",
        )
        job = SimpleNamespace(
            inputfile="/tmp/mol.xyz",
            pymol_script="/no/such/file.py",
            label="mol",
            quiet_mode=False,
            command_line_only=False,
        )
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        with pytest.raises(AssertionError, match="does not exist"):
            runner._get_visualization_command(job)

    def test_setup_style_default_uses_existing_cwd_style_file(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "zhang_group_pymol_style.py").write_text("# stub\n")
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(style=None, label="mol")
        command = runner._setup_style(job, "cmd")
        assert command == 'cmd -d "pymol_style mol'

    @pytest.mark.parametrize(
        "style,expected",
        [
            ("cylview", 'cmd -d "cylview_style mol'),
            ("cylview-flat", 'cmd -d "cylview_flat_style mol'),
        ],
    )
    def test_setup_style_cylview_variants(self, style, expected):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(style=style, label="mol")
        assert runner._setup_style(job, "cmd") == expected

    def test_setup_style_invalid_style_raises(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(style="bogus_style", label="mol")
        with pytest.raises(ValueError, match="not available"):
            runner._setup_style(job, "cmd")

    def test_add_vdw_appends_when_requested(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(vdw=True, label="mol")
        assert runner._add_vdw(job, "cmd") == "cmd; add_vdw mol"

    def test_add_vdw_noop_when_not_requested(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(vdw=False, label="mol")
        assert runner._add_vdw(job, "cmd") == "cmd"

    def test_add_coordinates_labels_handles_angles_and_dihedrals(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(coordinates=[[1, 2, 3], [1, 2, 3, 4]])
        command = runner._add_coordinates_labels(job, "cmd")
        assert "angle a1, id 1, id 2, id 3" in command
        assert "dihedral di1, id 1, id 2, id 3, id 4" in command

    def test_offset_labels_sets_position_when_given(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(label_offset="(0,0,1.2)")
        command = runner._offset_labels(job, "cmd")
        assert command == "cmd; set label_position, (0,0,1.2)"

    def test_offset_labels_noop_when_none(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(label_offset=None)
        assert runner._offset_labels(job, "cmd") == "cmd"

    def test_add_ray_command_appends_when_trace_true(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(trace=True)
        assert runner._add_ray_command(job, "cmd") == "cmd; ray 2400,1800"

    def test_add_ray_command_noop_when_trace_false(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        job = SimpleNamespace(trace=False)
        assert runner._add_ray_command(job, "cmd") == "cmd"

    def test_job_specific_commands_base_passthrough(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        assert runner._job_specific_commands(SimpleNamespace(), "cmd") == "cmd"

    def test_hide_labels_appends_command(self):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        assert (
            runner._hide_labels(SimpleNamespace(), "cmd") == "cmd; hide labels"
        )

    def test_create_process_raises_on_nonzero_returncode(
        self, tmp_path, mocker
    ):
        runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        runner.running_directory = str(tmp_path)
        job = SimpleNamespace(
            errfile=str(tmp_path / "j.err"), logfile=str(tmp_path / "j.log")
        )
        mock_process = mocker.MagicMock()
        mock_process.wait.return_value = 1
        mock_process.returncode = 1
        mocker.patch(
            "chemsmart.jobs.mol.runner.subprocess.Popen",
            return_value=mock_process,
        )

        with pytest.raises(subprocess.CalledProcessError):
            runner._create_process(job, "echo hi", {})

    def test_write_hybrid_pml_overwrites_existing_file(self, tmp_path, mocker):
        runner = PyMOLHybridVisualizationJobRunner.__new__(
            PyMOLHybridVisualizationJobRunner
        )
        job = SimpleNamespace(folder=str(tmp_path), label="mol")
        (tmp_path / "mol.pml").write_text("# existing\n")

        mocker.patch.object(
            PyMOLHybridVisualizationJobRunner, "_write_default_pymol_style"
        )
        mocker.patch.object(
            PyMOLHybridVisualizationJobRunner, "_write_faded_colors"
        )
        mocker.patch.object(
            PyMOLHybridVisualizationJobRunner, "_write_highlighted_colors"
        )
        mocker.patch.object(
            PyMOLHybridVisualizationJobRunner, "_write_surface_settings"
        )

        result = runner._write_hybrid_pml(job)

        assert result == str(tmp_path / "mol.pml")

    def test_write_highlighted_colors_reuses_schemes_when_more_groups(self):
        """More groups than default color schemes triggers the
        reuse-with-multiplier warning branch."""
        runner = PyMOLHybridVisualizationJobRunner.__new__(
            PyMOLHybridVisualizationJobRunner
        )
        job = SimpleNamespace(
            groups=[f"{i}-{i + 4}" for i in range(0, 55, 5)],  # 11 groups
            colors=[],
            stick_radius=None,
        )
        buf = io.StringIO()
        runner._write_highlighted_colors(job, buf)
        content = buf.getvalue()
        # 11 groups > 10 default color schemes: colors are reused via
        # the multiplier, so group11 wraps back around to the first
        # color scheme (cbap).
        assert "util.cbap group11" in content

    def test_write_highlighted_colors_uses_custom_stick_radius(self):
        runner = PyMOLHybridVisualizationJobRunner.__new__(
            PyMOLHybridVisualizationJobRunner
        )
        job = SimpleNamespace(groups=["1-5"], colors=[], stick_radius=0.5)
        buf = io.StringIO()
        runner._write_highlighted_colors(job, buf)
        assert "set stick_radius, 0.5," in buf.getvalue()

    def test_write_surface_settings_uses_job_overrides(self):
        runner = PyMOLHybridVisualizationJobRunner.__new__(
            PyMOLHybridVisualizationJobRunner
        )
        job = SimpleNamespace(surface_color="blue", surface_transparency="0.5")
        buf = io.StringIO()
        runner._write_surface_settings(job, buf)
        content = buf.getvalue()
        assert "set surface_color, blue, all" in content
        assert "set transparency, 0.5, all" in content

    def test_write_surface_settings_uses_defaults(self):
        runner = PyMOLHybridVisualizationJobRunner.__new__(
            PyMOLHybridVisualizationJobRunner
        )
        job = SimpleNamespace(surface_color=None, surface_transparency=None)
        buf = io.StringIO()
        runner._write_surface_settings(job, buf)
        content = buf.getvalue()
        assert "set surface_color, grey, all" in content
        assert "set transparency, 0.7, all" in content


class TestPyMOLScientificStyleVisualizationJobRunnerHelpers:
    def test_format_style_command_raises_for_non_scientific_style(self):
        from chemsmart.jobs.mol.runner import (
            PyMOLScientificStyleVisualizationJobRunner,
        )

        job = SimpleNamespace(style="pymol")
        with pytest.raises(ValueError, match="not available"):
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, "mol"
            )

    def test_generate_visualization_style_script_overwrites_existing(
        self, tmp_path
    ):
        from chemsmart.jobs.mol.runner import (
            PyMOLScientificStyleVisualizationJobRunner,
        )

        runner = PyMOLScientificStyleVisualizationJobRunner.__new__(
            PyMOLScientificStyleVisualizationJobRunner
        )
        job = SimpleNamespace(style="comic", folder=str(tmp_path))
        dest = tmp_path / "zhang_group_scientific_styles.py"
        dest.write_text("# stub\n")

        result = runner._generate_visualization_style_script(job)

        assert result == str(dest)
        assert dest.read_text() != "# stub\n"


class TestPyMOLMovieJobRunnerHelpers:
    def test_setup_style_uses_cwd_style_file_when_present(
        self, tmp_path, monkeypatch
    ):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        monkeypatch.chdir(tmp_path)
        (tmp_path / "zhang_group_pymol_style.py").write_text("# stub\n")
        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        job = SimpleNamespace(job_basename="mol")
        assert runner._setup_style(job, "cmd") == 'cmd -d "movie_style mol'

    def test_set_ray_trace_frames_appends_when_trace_true(self):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        job = SimpleNamespace(trace=True)
        assert (
            runner._set_ray_trace_frames(job, "cmd")
            == "cmd; set ray_trace_frames, 1; set ray_trace_mode, 1"
        )

    def test_create_movie_skips_when_mp4_exists_and_no_overwrite(
        self, tmp_path
    ):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        (tmp_path / "mol.mp4").write_text("stub")
        job = SimpleNamespace(
            folder=str(tmp_path),
            job_basename="mol",
            outputfile=str(tmp_path / "mol.pse"),
            overwrite=False,
        )
        assert runner._create_movie(job) is None

    def test_create_movie_raises_when_no_png_frames(self, tmp_path):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        job = SimpleNamespace(
            folder=str(tmp_path),
            job_basename="mol",
            outputfile=str(tmp_path / "mol.pse"),
            overwrite=False,
        )
        with pytest.raises(FileNotFoundError, match="No PNG frames found"):
            runner._create_movie(job)

    def test_create_movie_overwrites_existing_mp4_and_cleans_up_pngs(
        self, tmp_path, mocker
    ):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        (tmp_path / "mol.mp4").write_text("stub")
        (tmp_path / "mol_frame_0001.png").write_text("stub")
        job = SimpleNamespace(
            folder=str(tmp_path),
            job_basename="mol",
            outputfile=str(tmp_path / "mol.pse"),
            overwrite=True,
        )
        mock_run = mocker.patch("chemsmart.jobs.mol.runner.subprocess.run")

        runner._create_movie(job)

        mock_run.assert_called_once()
        assert not (tmp_path / "mol_frame_0001.png").exists()

    def test_create_movie_ffmpeg_called_process_error_propagates(
        self, tmp_path, mocker
    ):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        (tmp_path / "mol_frame_0001.png").write_text("stub")
        job = SimpleNamespace(
            folder=str(tmp_path),
            job_basename="mol",
            outputfile=str(tmp_path / "mol.pse"),
            overwrite=False,
        )
        mocker.patch(
            "chemsmart.jobs.mol.runner.subprocess.run",
            side_effect=subprocess.CalledProcessError(
                1, "ffmpeg", stderr="boom"
            ),
        )

        with pytest.raises(subprocess.CalledProcessError):
            runner._create_movie(job)

    def test_create_movie_ffmpeg_not_installed_raises(self, tmp_path, mocker):
        from chemsmart.jobs.mol.runner import PyMOLMovieJobRunner

        runner = PyMOLMovieJobRunner.__new__(PyMOLMovieJobRunner)
        (tmp_path / "mol_frame_0001.png").write_text("stub")
        job = SimpleNamespace(
            folder=str(tmp_path),
            job_basename="mol",
            outputfile=str(tmp_path / "mol.pse"),
            overwrite=False,
        )
        mocker.patch(
            "chemsmart.jobs.mol.runner.subprocess.run",
            side_effect=FileNotFoundError("no ffmpeg"),
        )

        with pytest.raises(FileNotFoundError, match="FFmpeg not found"):
            runner._create_movie(job)


class TestPyMOLIRCMovieJobRunnerHelpers:
    def test_get_rotation_command_is_noop(self):
        from chemsmart.jobs.mol.runner import PyMOLIRCMovieJobRunner

        runner = PyMOLIRCMovieJobRunner.__new__(PyMOLIRCMovieJobRunner)
        assert runner._get_rotation_command(SimpleNamespace(), "cmd") == "cmd"


class TestPyMOLNCIJobRunnerHelpers:
    @pytest.mark.parametrize(
        "binary,intermediate,expected_fragment",
        [
            (True, False, "nci_binary mol"),
            (False, True, "nci_intermediate mol"),
            (False, False, "nci mol"),
        ],
    )
    def test_run_nci_command_modes(
        self, binary, intermediate, expected_fragment
    ):
        from chemsmart.jobs.mol.runner import PyMOLNCIJobRunner

        runner = PyMOLNCIJobRunner.__new__(PyMOLNCIJobRunner)
        job = SimpleNamespace(
            binary=binary, intermediate=intermediate, source_basename="mol"
        )
        command = runner._run_nci_command(job, "cmd")
        assert command == f"cmd; {expected_fragment}"


class TestPyMOLMOJobRunnerHelpers:
    def test_generate_mo_cube_file_no_selection_raises(self, mocker):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        mocker.patch.object(
            PyMOLMOJobRunner,
            "_get_gaussian_executable",
            return_value="/opt/g16",
        )
        job = SimpleNamespace(
            number=None,
            homo=False,
            lumo=False,
            job_basename="mol",
            source_basename="mol",
        )
        with pytest.raises(ValueError, match="exactly one of"):
            runner._generate_mo_cube_file(job)

    def test_generate_mo_cube_file_multiple_selections_raises(self, mocker):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        mocker.patch.object(
            PyMOLMOJobRunner,
            "_get_gaussian_executable",
            return_value="/opt/g16",
        )
        job = SimpleNamespace(
            number=5,
            homo=True,
            lumo=False,
            job_basename="mol",
            source_basename="mol",
        )
        with pytest.raises(ValueError, match="exactly one of"):
            runner._generate_mo_cube_file(job)

    def test_generate_mo_cube_file_skips_when_cube_exists(
        self, tmp_path, monkeypatch, mocker
    ):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        monkeypatch.chdir(tmp_path)
        (tmp_path / "mol.cube").write_text("stub")
        mocker.patch.object(
            PyMOLMOJobRunner,
            "_get_gaussian_executable",
            return_value="/opt/g16",
        )
        mock_run = mocker.patch("chemsmart.jobs.mol.runner.run_command")
        job = SimpleNamespace(
            number=None,
            homo=True,
            lumo=False,
            job_basename="mol",
            source_basename="mol",
        )

        assert runner._generate_mo_cube_file(job) is None
        mock_run.assert_not_called()

    @pytest.mark.parametrize(
        "kwargs,expected_mo",
        [
            ({"number": 7, "homo": False, "lumo": False}, "7"),
            ({"number": None, "homo": True, "lumo": False}, "HOMO"),
            ({"number": None, "homo": False, "lumo": True}, "LUMO"),
        ],
    )
    def test_generate_mo_cube_file_runs_cubegen(
        self, tmp_path, monkeypatch, mocker, kwargs, expected_mo
    ):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        monkeypatch.chdir(tmp_path)
        mocker.patch.object(
            PyMOLMOJobRunner,
            "_get_gaussian_executable",
            return_value="/opt/g16",
        )
        mock_run = mocker.patch("chemsmart.jobs.mol.runner.run_command")
        job = SimpleNamespace(
            job_basename="mol", source_basename="mol", **kwargs
        )

        runner._generate_mo_cube_file(job)

        assert mock_run.call_args.args[0] == (
            f"/opt/g16/cubegen 0 MO={expected_mo} mol.fchk mol.cube 0 h"
        )

    def test_write_molecular_orbital_pml_overwrites_existing(self, tmp_path):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        job = SimpleNamespace(
            folder=str(tmp_path),
            mo_basename="mol_HOMO",
            isosurface_value=0.05,
            transparency_value=0.3,
            surface_quality=1,
            antialias_value=2,
        )
        pml_path = tmp_path / "mol_HOMO.pml"
        pml_path.write_text("# existing\n")

        runner._write_molecular_orbital_pml(job)

        content = pml_path.read_text()
        assert "load mol_HOMO.cube" in content
        assert "isosurface pos_iso, mol_HOMO, 0.05" in content
        assert "isosurface neg_iso, mol_HOMO, -0.05" in content

    def test_offset_labels_is_noop(self):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        assert runner._offset_labels(SimpleNamespace(), "cmd") == "cmd"

    def test_call_pml_appends_load_command(self):
        from chemsmart.jobs.mol.runner import PyMOLMOJobRunner

        runner = PyMOLMOJobRunner.__new__(PyMOLMOJobRunner)
        job = SimpleNamespace(folder="/tmp/testjob", mo_basename="mol_HOMO")
        command = runner._call_pml(job, "cmd")
        assert command == "cmd; load /tmp/testjob/mol_HOMO.pml"


class TestPyMOLSpinJobRunnerHelpers:
    def test_get_gaussian_executable(self, mocker):
        from chemsmart.jobs.mol.runner import PyMOLSpinJobRunner

        runner = PyMOLSpinJobRunner.__new__(PyMOLSpinJobRunner)
        runner.server = SimpleNamespace(name="local")
        mock_exe_cls = mocker.patch(
            "chemsmart.jobs.mol.runner.GaussianExecutable"
        )
        mock_exe = mocker.MagicMock()
        mock_exe.executable_folder = "/opt/g16"
        mock_exe_cls.from_servername.return_value = mock_exe

        result = runner._get_gaussian_executable(SimpleNamespace())

        assert result == "/opt/g16"
        mock_exe_cls.from_servername.assert_called_once_with("local")

    def test_generate_spin_cube_file_runs_cubegen(self, mocker):
        from chemsmart.jobs.mol.runner import PyMOLSpinJobRunner

        runner = PyMOLSpinJobRunner.__new__(PyMOLSpinJobRunner)
        mocker.patch.object(
            PyMOLSpinJobRunner,
            "_get_gaussian_executable",
            return_value="/opt/g16",
        )
        mock_run = mocker.patch("chemsmart.jobs.mol.runner.run_command")
        job = SimpleNamespace(
            source_basename="mol", job_basename="mol", npts=100
        )

        runner._generate_spin_cube_file(job)

        assert (
            mock_run.call_args.args[0]
            == "/opt/g16/cubegen 0 spin mol.fchk mol.cube 100"
        )

    def test_write_spin_density_pml_overwrites_existing(self, tmp_path):
        from chemsmart.jobs.mol.runner import PyMOLSpinJobRunner

        runner = PyMOLSpinJobRunner.__new__(PyMOLSpinJobRunner)
        job = SimpleNamespace(
            folder=str(tmp_path),
            spin_basename="mol_spin",
            isosurface_value=0.004,
            transparency_value=0.3,
            surface_quality=1,
            antialias_value=2,
            ray_trace_mode=1,
        )
        pml_path = tmp_path / "mol_spin.pml"
        pml_path.write_text("# existing\n")

        runner._write_spin_density_pml(job)

        content = pml_path.read_text()
        assert "load mol_spin.cube" in content
        assert "isosurface pos_iso_spin, mol_spin, 0.004" in content
        assert "isosurface neg_iso_spin, mol_spin, -0.004" in content
        assert "set ray_trace_mode, 1" in content

    def test_job_specific_commands_chains_hide_pml_and_ray(self, mocker):
        from chemsmart.jobs.mol.runner import PyMOLSpinJobRunner

        runner = PyMOLSpinJobRunner.__new__(PyMOLSpinJobRunner)
        job = SimpleNamespace(
            folder="/tmp/testjob", spin_basename="mol_spin", trace=True
        )
        command = runner._job_specific_commands(job, "cmd")
        assert "hide labels" in command
        assert "load /tmp/testjob/mol_spin.pml" in command
        assert "ray 2400,1800" in command

    def test_offset_labels_is_noop(self):
        from chemsmart.jobs.mol.runner import PyMOLSpinJobRunner

        runner = PyMOLSpinJobRunner.__new__(PyMOLSpinJobRunner)
        assert runner._offset_labels(SimpleNamespace(), "cmd") == "cmd"

    def test_call_pml_appends_load_command(self):
        from chemsmart.jobs.mol.runner import PyMOLSpinJobRunner

        runner = PyMOLSpinJobRunner.__new__(PyMOLSpinJobRunner)
        job = SimpleNamespace(folder="/tmp/testjob", spin_basename="mol_spin")
        command = runner._call_pml(job, "cmd")
        assert command == "cmd; load /tmp/testjob/mol_spin.pml"


class TestPyMOLAlignJobRunnerWriteInputAndRun:
    """Direct tests for PyMOLAlignJobRunner._write_input and run."""

    def _make_molecule(self, name, x=0.0):
        mol = Molecule(
            symbols=["Ar"], positions=[[x, 0.0, 0.0]], charge=0, multiplicity=1
        )
        mol.name = name
        return mol

    def test_write_input_writes_xyz_files_and_sets_no_batch(self, tmp_path):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        mols = [self._make_molecule("mol1"), self._make_molecule("mol2", 1.0)]
        job = SimpleNamespace(folder=str(tmp_path), molecule=mols)

        runner._write_input(job)

        assert job.use_batch_processing is False
        assert job.total_batches == 1
        assert job.mol_names == ["mol1", "mol2"]
        assert all(os.path.exists(p) for p in job.xyz_absolute_paths)

    def test_write_input_skips_writing_when_xyz_already_exists(self, tmp_path):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        (tmp_path / "mol1.xyz").write_text("existing stub\n")
        job = SimpleNamespace(
            folder=str(tmp_path), molecule=[self._make_molecule("mol1")]
        )

        runner._write_input(job)

        # file content is untouched (not overwritten by mol.write)
        assert (tmp_path / "mol1.xyz").read_text() == "existing stub\n"

    def test_write_input_non_molecule_raises(self, tmp_path):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        job = SimpleNamespace(
            folder=str(tmp_path), molecule=["not-a-molecule"]
        )
        with pytest.raises(ValueError, match="not of Molecule type"):
            runner._write_input(job)

    def test_write_input_missing_name_attribute_raises(self, tmp_path):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        mol = Molecule(
            symbols=["Ar"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        job = SimpleNamespace(folder=str(tmp_path), molecule=[mol])
        with pytest.raises(ValueError, match="missing .name attribute"):
            runner._write_input(job)

    def test_write_input_enables_batch_processing_above_threshold(
        self, tmp_path
    ):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        runner.MAX_MOLECULES_PER_BATCH = 3
        mols = [self._make_molecule(f"m{i}", float(i)) for i in range(7)]
        job = SimpleNamespace(folder=str(tmp_path), molecule=mols)

        runner._write_input(job)

        assert job.use_batch_processing is True
        assert job.total_batches == 3

    def test_run_dispatches_to_batch_processing(self, mocker):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)

        def fake_write_input(job):
            job.use_batch_processing = True
            job.total_batches = 2
            job.xyz_absolute_paths = ["a.xyz", "b.xyz"]

        mocker.patch.object(PyMOLAlignJobRunner, "_prerun")
        mocker.patch.object(
            PyMOLAlignJobRunner, "_write_input", side_effect=fake_write_input
        )
        mock_batch = mocker.patch.object(
            PyMOLAlignJobRunner,
            "_run_batch_processing",
            return_value="batch-result",
        )

        job = SimpleNamespace()
        result = runner.run(job)

        assert result == "batch-result"
        mock_batch.assert_called_once_with(job)

    def test_run_dispatches_to_single_batch_path(self, mocker):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)

        def fake_write_input(job):
            job.use_batch_processing = False
            job.xyz_absolute_paths = ["a.xyz"]

        mocker.patch.object(PyMOLAlignJobRunner, "_prerun")
        mocker.patch.object(
            PyMOLAlignJobRunner, "_write_input", side_effect=fake_write_input
        )
        mocker.patch.object(
            PyMOLAlignJobRunner, "_get_command", return_value="cmd"
        )
        mocker.patch.object(
            PyMOLAlignJobRunner, "_update_os_environ", return_value={}
        )
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")
        mock_postrun = mocker.patch.object(PyMOLAlignJobRunner, "_postrun")
        mock_cleanup = mocker.patch.object(
            PyMOLAlignJobRunner, "_postrun_cleanup"
        )

        job = SimpleNamespace()
        result = runner.run(job)

        assert result is job
        mock_postrun.assert_called_once_with(job)
        mock_cleanup.assert_called_once_with(job)


class TestPyMOLAlignJobRunnerBatchProcessing:
    """Direct tests for PyMOLAlignJobRunner._run_batch_processing and
    _execute_batch, exercised via a bare instance with the subprocess-
    launching methods (_create_process/_run) mocked out, since these
    build pure command strings and never touch pymol.cmd directly."""

    def _make_runner(self):
        from unittest.mock import PropertyMock, patch

        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        runner.running_directory = "/tmp"
        self._executable_patcher = patch.object(
            PyMOLAlignJobRunner,
            "executable",
            new_callable=PropertyMock,
            return_value="/usr/bin/pymol",
        )
        self._executable_patcher.start()
        return runner

    def _make_job(self, n_molecules=3, style=None, total_batches=1):
        return SimpleNamespace(
            folder="/tmp/testjob",
            label="testalign",
            total_batches=total_batches,
            xyz_absolute_paths=[f"/tmp/m{i}.xyz" for i in range(n_molecules)],
            mol_names=[f"mol{i}" for i in range(n_molecules)],
            style=style,
            quiet_mode=False,
            command_line_only=False,
            pymol_script="/tmp/style.py",
        )

    def test_run_batch_processing_executes_one_batch_per_chunk(self, mocker):
        runner = self._make_runner()
        job = self._make_job(n_molecules=250, total_batches=3)
        runner.MAX_MOLECULES_PER_BATCH = 100

        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mock_create = mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(returncode=0),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")
        mock_postrun = mocker.patch.object(PyMOLAlignJobRunner, "_postrun")
        mock_cleanup = mocker.patch.object(
            PyMOLAlignJobRunner, "_postrun_cleanup"
        )

        result = runner._run_batch_processing(job)

        assert result is job
        assert mock_create.call_count == 3
        mock_postrun.assert_called_once_with(job)
        mock_cleanup.assert_called_once_with(job)
        # First batch loads molecules directly (no append mode); later
        # batches append to the same log/err files.
        append_modes = [
            c.kwargs.get("append_mode") for c in mock_create.call_args_list
        ]
        assert append_modes == [False, True, True]

    def test_execute_batch_first_batch_loads_molecules_directly(self, mocker):
        runner = self._make_runner()
        job = self._make_job(n_molecules=2)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mock_create = mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(returncode=0),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")

        runner._execute_batch(job, 0, job.xyz_absolute_paths, job.mol_names)

        command = mock_create.call_args.args[1]
        assert command.startswith("/usr/bin/pymol /tmp/m0.xyz /tmp/m1.xyz")
        assert "pymol_style mol0" in command
        assert "pymol_style mol1" in command
        assert "align mol1, mol0" in command
        assert f"save {job.folder}/testalign.pse".replace(
            "/tmp/testjob", job.folder
        )
        assert mock_create.call_args.kwargs["append_mode"] is False

    def test_execute_batch_subsequent_batch_opens_existing_pse(self, mocker):
        runner = self._make_runner()
        job = self._make_job(n_molecules=3)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mock_create = mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(returncode=0),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")

        runner._execute_batch(
            job, 1, job.xyz_absolute_paths[1:], ["mol1", "mol2"]
        )

        command = mock_create.call_args.args[1]
        assert command.startswith("/usr/bin/pymol /tmp/testjob/testalign.pse")
        assert "load /tmp/m1.xyz" in command
        assert "load /tmp/m2.xyz" in command
        # subsequent batches align every molecule to the global first
        # molecule, without excluding it by name (unlike batch 0)
        assert "align mol1, mol0" in command
        assert "align mol2, mol0" in command
        assert mock_create.call_args.kwargs["append_mode"] is True

    @pytest.mark.parametrize(
        "style,expected_cmd",
        [
            (None, "pymol_style"),
            ("pymol", "pymol_style"),
            ("cylview", "cylview_style"),
            ("cylview-flat", "cylview_flat_style"),
        ],
    )
    def test_execute_batch_style_variants(self, mocker, style, expected_cmd):
        runner = self._make_runner()
        job = self._make_job(n_molecules=2, style=style)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mock_create = mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(returncode=0),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")

        runner._execute_batch(job, 0, job.xyz_absolute_paths, job.mol_names)

        command = mock_create.call_args.args[1]
        assert f"{expected_cmd} mol0" in command
        assert f"{expected_cmd} mol1" in command

    def test_execute_batch_invalid_style_raises(self, mocker):
        runner = self._make_runner()
        job = self._make_job(n_molecules=2, style="bogus_style")
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )

        with pytest.raises(ValueError, match="not available"):
            runner._execute_batch(
                job, 0, job.xyz_absolute_paths, job.mol_names
            )

    def test_execute_batch_nonzero_returncode_raises_runtime_error(
        self, mocker
    ):
        runner = self._make_runner()
        job = self._make_job(n_molecules=2)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(returncode=1),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")

        with pytest.raises(RuntimeError, match="execution failed"):
            runner._execute_batch(
                job, 0, job.xyz_absolute_paths, job.mol_names
            )

    def test_execute_batch_propagates_and_logs_unexpected_exception(
        self, mocker
    ):
        runner = self._make_runner()
        job = self._make_job(n_molecules=2)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            side_effect=OSError("simulated launch failure"),
        )

        with pytest.raises(OSError, match="simulated launch failure"):
            runner._execute_batch(
                job, 0, job.xyz_absolute_paths, job.mol_names
            )

    def test_quiet_mode_and_command_line_only_flags_appended(self, mocker):
        runner = self._make_runner()
        job = self._make_job(n_molecules=1)
        job.quiet_mode = True
        job.command_line_only = True
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        mock_create = mocker.patch.object(
            PyMOLAlignJobRunner,
            "_create_process",
            return_value=mocker.MagicMock(returncode=0),
        )
        mocker.patch.object(PyMOLAlignJobRunner, "_run")

        runner._execute_batch(job, 0, job.xyz_absolute_paths, job.mol_names)

        command = mock_create.call_args.args[1]
        assert " -q" in command
        assert " -c" in command


class TestPyMOLAlignJobRunnerNonBatchCommandHelpers:
    """Direct tests for PyMOLAlignJobRunner's own overrides of
    _add_style_script, _get_visualization_command, _setup_style, and
    _align_command (used by the non-batch-processing path)."""

    def test_add_style_script_uses_existing_folder_style_file(self, tmp_path):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        style_file = tmp_path / "zhang_group_pymol_style.py"
        style_file.write_text("# stub\n")
        job = SimpleNamespace(pymol_script=None, folder=str(tmp_path))

        command = runner._add_style_script(job, "cmd")

        assert command == f"cmd -r {style_file}"

    def test_add_style_script_uses_existing_user_script(self, tmp_path):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        script = tmp_path / "custom_style.py"
        script.write_text("# stub\n")
        job = SimpleNamespace(pymol_script=str(script))

        command = runner._add_style_script(job, "cmd")

        assert command == f"cmd -r {script}"

    def test_add_style_script_missing_user_script_raises(self):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        job = SimpleNamespace(pymol_script="/no/such/file.py")
        with pytest.raises(FileNotFoundError, match="does not exist"):
            runner._add_style_script(job, "cmd")

    def test_get_visualization_command_raises_when_no_xyz_files(self, mocker):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "executable",
            new_callable=mocker.PropertyMock,
            return_value="/usr/bin/pymol",
        )
        job = SimpleNamespace(xyz_absolute_paths=[])
        with pytest.raises(ValueError, match="No XYZ files found"):
            runner._get_visualization_command(job)

    def test_get_visualization_command_single_file_with_flags(self, mocker):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "executable",
            new_callable=mocker.PropertyMock,
            return_value="/usr/bin/pymol",
        )
        mocker.patch.object(
            PyMOLAlignJobRunner,
            "_add_style_script",
            side_effect=lambda job, cmd: cmd,
        )
        job = SimpleNamespace(
            xyz_absolute_paths=["/tmp/m0.xyz"],
            quiet_mode=True,
            command_line_only=True,
        )

        command = runner._get_visualization_command(job)

        assert command == "/usr/bin/pymol /tmp/m0.xyz -q -c"

    def test_setup_style_cylview(self):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        job = SimpleNamespace(style="cylview", mol_names=["mol1", "mol2"])
        command = runner._setup_style(job, "cmd")
        assert command == 'cmd -d "cylview_style mol1; cylview_style mol2'

    def test_setup_style_invalid_raises(self):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        job = SimpleNamespace(style="bogus_style", mol_names=["mol1"])
        with pytest.raises(ValueError, match="not available"):
            runner._setup_style(job, "cmd")

    def test_align_command_noop_with_single_molecule(self):
        runner = PyMOLAlignJobRunner.__new__(PyMOLAlignJobRunner)
        job = SimpleNamespace(mol_names=["mol1"])
        assert runner._align_command(job, "cmd") == "cmd"


@pytest.mark.usefixtures("skip_if_no_pymol")
class TestPyMOLStyleCommandsContinued:
    label_1_mer = "1-mer"

    def test_scientific_style_registry_matches_classes(self):
        template_commands = (
            zhang_group_scientific_styles.PYMOL_SCIENTIFIC_STYLE_COMMANDS
        )

        assert PYMOL_SCIENTIFIC_STYLE_COMMANDS == template_commands
        assert len(SCIENTIFIC_STYLE_CLASSES) == len(template_commands)
        assert PYMOL_VISUALIZE_STYLE_CLI_CHOICES == [
            style_cls.command.replace("_", "-")
            for style_cls in SCIENTIFIC_STYLE_CLASSES
        ]
        for style_cls in SCIENTIFIC_STYLE_CLASSES:
            assert template_commands[style_cls.command] == style_cls.command
            wrapper = getattr(zhang_group_scientific_styles, style_cls.command)
            assert wrapper.__name__ == style_cls.command

    def test_select_coordination_uses_command_as_prefix(self, mocker):
        style = StericSurfaceStyle()
        mock_build = mocker.patch.object(
            ScientificStyle, "build_coordination_atoms", return_value={}
        )
        style.select_coordination("all")
        mock_build.assert_called_once_with(
            selection="all",
            prefix="steric_surface",
            include_nh_h=False,
        )

    def test_comic_metallic_reuses_glossy_palette(self):
        assert ComicMetallicStyle.colors is GlossyStyle.colors
        assert ComicMetallicStyle.DONOR_ROLES

    def test_zhang_group_scientific_styles_hide_distance_value_labels(
        self, mocker
    ):
        mock_cmd = mocker.patch(
            "chemsmart.jobs.mol.templates.zhang_group_scientific_styles.cmd"
        )
        mocker.patch.object(
            ScientificStyle,
            "_distance_object_names",
            return_value=["d1", "d2"],
        )
        ScientificStyle.hide_distance_value_labels()
        mock_cmd.hide.assert_any_call("labels", "d1")
        mock_cmd.hide.assert_any_call("labels", "d2")
        assert mock_cmd.hide.call_count == 2

    def test_finish_default_runs_finalize_before_camera(self, mocker):
        style = ComicMetallicStyle()
        calls = []

        def record_finalize():
            calls.append("finalize")

        def record_camera(selection):
            calls.append(("camera", selection))

        mocker.patch.object(style, "finalize", side_effect=record_finalize)
        mocker.patch.object(style, "finish_camera", side_effect=record_camera)
        style.finish_default("all")
        assert calls == ["finalize", ("camera", "all")]

    def test_scientific_style_runner_orders_distances_before_render(self):
        runner = PyMOLScientificStyleVisualizationJobRunner.__new__(
            PyMOLScientificStyleVisualizationJobRunner
        )
        job = SimpleNamespace(
            style="comic",
            coordinates=[[1, 8]],
            label=self.label_1_mer,
        )
        command = runner._setup_style(job, "pymol cmd")
        command = runner._add_coordinates_labels(job, command)
        command = runner._append_style_render(job, command)
        assert command.index("distance d1, id 1, id 8") < command.index(
            "comic"
        )
        assert "hide labels" not in command

    def test_style_wrapper_ignores_pymol_self_kwarg(self, mocker):
        mock_render = mocker.patch.object(ComicMetallicStyle, "render")
        getattr(zhang_group_scientific_styles, "comic")("sel", _self="ignored")
        mock_render.assert_called_once_with(selection="sel")

    def test_build_coordination_atoms_uses_geometry_helpers(self):
        assert hasattr(zhang_group_scientific_styles, "get_coordinating_atoms")
        assert hasattr(zhang_group_scientific_styles, "is_metal")
        assert hasattr(ScientificStyle, "_role_pymol_indices")

    def test_element_category_helpers_are_centralized(self):
        assert ScientificStyle.ELEMENT_CATEGORIES["halogen"] == "F+Cl+Br+I"
        assert ScientificStyle.ELEMENT_CATEGORIES["N+O"] == "N+O"
        assert ScientificStyle.ELEMENT_CATEGORIES["metal"].startswith("elem ")
        assert (
            ScientificStyle.element_category_selection("all", "halogen")
            == "(all) and elem F+Cl+Br+I"
        )
        assert (
            ScientificStyle.element_category_selection("sc_shell", "N+O")
            == "(sc_shell) and elem N+O"
        )
        assert hasattr(ScientificStyle, "apply_element_palette")
        assert hasattr(ScientificStyle, "_safe_set")

    def test_pymol_elem_selection(self):
        assert pymol_elem_selection(["Fe", "Cu"]) == "elem Fe+Cu"
        assert pymol_elem_selection([]) == "none"

    def test_metal_pymol_selection(self):
        selection = metal_pymol_selection()
        assert selection.startswith("elem ")
        symbols = selection.replace("elem ", "").split("+")
        assert "Mn" in symbols
        assert "C" not in symbols
        assert selection == ScientificStyle.METAL_ELEMENTS

    def test_hybrid_is_not_a_derived_style(self):
        # TODO: Remove this test if hybrid style is integrated with ScientificStyle
        with pytest.raises(ValueError, match="not available"):
            normalize_pymol_style("hybrid")

    def test_derived_styles_hide_distance_value_labels_base_style_keeps_them(
        self,
    ):
        base_runner = PyMOLJobRunner.__new__(PyMOLJobRunner)
        derived_runner = PyMOLScientificStyleVisualizationJobRunner.__new__(
            PyMOLScientificStyleVisualizationJobRunner
        )
        job = SimpleNamespace(style="comic", coordinates=[[1, 8]])

        base_command = base_runner._add_coordinates_labels(job, "cmd")
        derived_command = derived_runner._add_coordinates_labels(job, "cmd")

        assert "distance d1, id 1, id 8" in base_command
        assert "distance d1, id 1, id 8" in derived_command
        assert "hide labels" not in base_command
        assert "hide labels" not in derived_command
        assert hasattr(ScientificStyle, "hide_distance_value_labels")


@pytest.mark.usefixtures("skip_if_no_pymol")
class TestPyMOLScientificStyleVisualizationJobs:
    """Derived ``-s`` style jobs using the ``1-mer.xyz`` test geometry."""

    coordination_bonds_1_mer = [
        [1, 2],
        [1, 5],
        [1, 36],
        [1, 3],
        [1, 15],
        [1, 8],
    ]

    def test_comic_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="comic",
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"comic {job.label}"
        )

    def test_comic_style_job_with_coordination_bonds_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="comic",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"comic {job.label}"
        )
        runner = PyMOLScientificStyleVisualizationJobRunner.__new__(
            PyMOLScientificStyleVisualizationJobRunner
        )
        assert "distance d1, id 1, id 2" in runner._add_coordinates_labels(
            job, "cmd"
        )
        assert "hide labels" not in runner._add_coordinates_labels(job, "cmd")

    def test_glossy_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="glossy",
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"glossy {job.label}"
        )

    def test_editorial_minimal_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="editorial-minimal",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"editorial_minimal {job.label}"
        )

    def test_soft_ceramic_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="soft-ceramic",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"soft_ceramic {job.label}"
        )

    def test_soft_cartoon_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="soft-cartoon",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"soft_cartoon {job.label}"
        )

    def test_neon_coordination_core_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="neon-coordination-core",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"neon_coordination_core {job.label}"
        )

    def test_matte_clay_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="matte-clay",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"matte_clay {job.label}"
        )

    def test_xray_wire_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="xray-wire",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"xray_wire {job.label}"
        )

    def test_steric_surface_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="steric-surface",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"steric_surface {job.label}"
        )

    def test_quasi_chemdraw_bold_style_job_on_1_mer_xyz(
        self,
        tmpdir,
        visualized_1_mer_xyz_file,
        pymol_scientific_style_visualization_jobrunner,
    ):
        job = PyMOLScientificStyleVisualizationJob.from_filename(
            visualized_1_mer_xyz_file,
            jobrunner=pymol_scientific_style_visualization_jobrunner,
            style="quasi-chemdraw-bold",
            coordinates=self.coordination_bonds_1_mer,
        )
        job.set_folder(tmpdir)
        job.run()

        assert job.is_complete()
        assert os.path.exists(
            os.path.join(tmpdir, "zhang_group_scientific_styles.py")
        )
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"quasi_chemdraw_bold {job.label}"
        )


@pytest.fixture()
def loaded_1_mer(visualized_1_mer_xyz_file):
    """Load ``1-mer.xyz`` into a fresh, real PyMOL session (PyMOL required)."""
    from pymol import cmd

    cmd.reinitialize()
    cmd.load(visualized_1_mer_xyz_file, "mol1mer")
    yield "mol1mer"
    cmd.delete("all")


@pytest.mark.usefixtures("skip_if_no_pymol")
class TestPyMOLScientificStyleRenderDirect:
    """Call ``ScientificStyle.render()`` in-process against a real PyMOL
    session so ``coverage`` observes the style bodies. ``job.run()`` (used
    in :class:`TestPyMOLScientificStyleVisualizationJobs`) shells out to a
    ``pymol`` subprocess, so it never exercises this module's lines in the
    test process itself.
    """

    @pytest.mark.parametrize(
        "style_cls",
        SCIENTIFIC_STYLE_CLASSES,
        ids=[style_cls.command for style_cls in SCIENTIFIC_STYLE_CLASSES],
    )
    def test_render_each_style_directly_on_real_pymol_session(
        self, style_cls, loaded_1_mer
    ):
        from pymol import cmd

        style_cls().render(loaded_1_mer)
        assert cmd.count_atoms(loaded_1_mer) > 0

    def test_comic_metallic_render_bonds_highlighted_distance_pairs(
        self, loaded_1_mer
    ):
        """``-c`` distance objects (``d1``, ``d2``, …) drive explicit
        metal-donor bonding in :meth:`ComicMetallicStyle.render`."""
        from pymol import cmd

        cmd.distance(
            "d1", f"{loaded_1_mer} and id 1", f"{loaded_1_mer} and id 2"
        )
        cmd.distance(
            "d2", f"{loaded_1_mer} and id 1", f"{loaded_1_mer} and id 5"
        )
        assert ScientificStyle.pairs_from_distance_objects(loaded_1_mer)

        ComicMetallicStyle().render(loaded_1_mer)

        # `render()` consumes and removes the highlight distance objects.
        assert ScientificStyle._distance_object_names() == []

    def test_select_coordination_handles_empty_selection(self, loaded_1_mer):
        """An empty selection short-circuits ``build_coordination_atoms``."""
        atoms = ScientificStyle.build_coordination_atoms(
            selection="none", prefix="empty_test"
        )
        from pymol import cmd

        for name in atoms.values():
            assert cmd.count_atoms(name) == 0

    def test_build_coordination_atoms_accepts_explicit_metal(
        self, loaded_1_mer
    ):
        """Passing ``metal=`` bypasses automatic ``is_metal`` detection."""
        from pymol import cmd

        atoms = ScientificStyle.build_coordination_atoms(
            selection=loaded_1_mer,
            prefix="explicit_metal_test",
            metal="elem Mn",
        )
        assert cmd.count_atoms(atoms["metal"]) > 0

    def test_apply_coordination_sci_palette_runs_on_real_selection(
        self, loaded_1_mer
    ):
        style = ScientificStyle()
        style.define_shared_colors()
        atoms = style.select_coordination(loaded_1_mer)
        style.apply_coordination_sci_palette(loaded_1_mer, atoms)

    def test_frame_and_finalize_on_real_selection(self, loaded_1_mer):
        style = ScientificStyle()
        atoms = style.select_coordination(loaded_1_mer)
        style.frame(loaded_1_mer, atoms["coordination_core"])
        style.finalize()
        style.finish_default(loaded_1_mer)

    def test_apply_camera_and_lighting_helpers_on_real_session(
        self, loaded_1_mer
    ):
        style = ScientificStyle()
        style.set_transparent_background()
        style.apply_base_quality()
        style.apply_lighting(
            0.5, 0.5, 0.3, 0.7, 0.2, spec_power=100, shininess=50
        )
        style.apply_transparent_view(
            orthoscopic=1,
            field_of_view=30,
            depth_cue=1,
            fog_start=0.5,
            ray_trace_gain=0.1,
        )
        style.apply_illustrated_camera()
        style.apply_soft_shadows()
        style.apply_ambient_occlusion()
        style.define_shared_colors()
        style.finish_camera(loaded_1_mer)
        style.safe_ray_shadows("light")

        # Optional-argument branches not exercised by any style's render():
        # `apply_lighting` without spec_power, and `apply_transparent_view`
        # without depth_cue / fog_start / ray_shadows_mode.
        style.apply_lighting(0.4, 0.2, 0.3, 0.6, 0.1)
        style.apply_transparent_view(
            orthoscopic=1,
            field_of_view=30,
            depth_cue=None,
            fog_start=None,
            ray_trace_gain=None,
            ray_shadows_mode=None,
        )

    def test_safe_set_and_pairs_helpers_with_no_distance_objects(
        self, loaded_1_mer
    ):
        style = ScientificStyle()
        style.safe_set("orthoscopic", 1)
        style._safe_set(
            "sphere_scale", 0.3, selection=loaded_1_mer, category="C"
        )
        # No selection/category: takes the direct `safe_set` branch.
        style._safe_set("orthoscopic", 1)
        assert ScientificStyle.pairs_from_distance_objects(loaded_1_mer) == []
        ScientificStyle.remove_distance_objects()
        ScientificStyle.hide_distance_value_labels()

    def test_coordination_sphere_atoms_without_hydride(self, loaded_1_mer):
        style = ScientificStyle()
        atoms = style.select_coordination(loaded_1_mer)
        with_hydride = style.coordination_sphere_atoms(atoms)
        without_hydride = style.coordination_sphere_atoms(
            atoms, include_hydride=False
        )
        assert isinstance(with_hydride, str)
        assert isinstance(without_hydride, str)

    def test_base_render_raises_not_implemented(self):
        with pytest.raises(NotImplementedError, match="must implement render"):
            ScientificStyle().render()

    def test_metal_element_label_returns_placeholder_for_empty_selection(
        self, loaded_1_mer
    ):
        assert GlossyStyle.metal_element_label("none") == "?"

    def test_element_category_selection_metal_uses_elem_prefixed_branch(self):
        selection = ScientificStyle.element_category_selection("all", "metal")
        assert selection.startswith("(all) and (elem ")

    def test_comic_metallic_render_without_metal_present(self, loaded_1_mer):
        """No metal in the rendered selection exercises the ``count_atoms
        == 0`` branches guarding bonding/labeling in
        :meth:`ComicMetallicStyle.render`."""
        from pymol import cmd

        non_metal_selection = f"{loaded_1_mer} and elem C+H"
        assert cmd.count_atoms(non_metal_selection) > 0
        ComicMetallicStyle().render(non_metal_selection)

    def test_steric_surface_render_skips_spheres_when_none_are_selected(
        self, loaded_1_mer, mocker
    ):
        """``select_coordination`` always returns named PyMOL selections
        (never a literal ``"none"`` string), so ``coordination_sphere_atoms``
        is always truthy through a real render. Force it empty directly to
        exercise the ``if sphere_atoms:`` False branch in
        :meth:`StericSurfaceStyle.render`."""
        mocker.patch.object(
            StericSurfaceStyle, "coordination_sphere_atoms", return_value=""
        )
        StericSurfaceStyle().render(loaded_1_mer)

    def test_parse_metal_symbols_returns_empty_set_for_falsy_input(self):
        assert ScientificStyle._parse_metal_symbols(None) == set()
        assert ScientificStyle._parse_metal_symbols("") == set()

    def test_build_coordination_atoms_lone_metal_has_no_donors(self, tmp_path):
        """A metal with no nearby atoms exercises the ``donors -> none``
        branch (no donor_s/donor_n/donor_p atoms found)."""
        from pymol import cmd

        xyz_path = tmp_path / "lone_metal.xyz"
        xyz_path.write_text("1\nlone metal\nMn 0.0 0.0 0.0\n")
        cmd.reinitialize()
        cmd.load(str(xyz_path), "lone_metal")
        try:
            atoms = ScientificStyle.build_coordination_atoms(
                selection="lone_metal", prefix="lone_metal_test"
            )
            assert cmd.count_atoms(atoms["metal"]) == 1
            assert cmd.count_atoms(atoms["donors"]) == 0
        finally:
            cmd.delete("all")

    def test_build_coordination_atoms_donor_p_primary_shell(self, tmp_path):
        """A phosphorus donor within primary bonding range of the metal
        exercises the ``donor_p_pymol`` branch in ``donor_parts``."""
        from pymol import cmd

        xyz_path = tmp_path / "metal_phosphine.xyz"
        xyz_path.write_text(
            "2\nmetal + P donor\nMn 0.0 0.0 0.0\nP 0.0 0.0 2.2\n"
        )
        cmd.reinitialize()
        cmd.load(str(xyz_path), "metal_p")
        try:
            atoms = ScientificStyle.build_coordination_atoms(
                selection="metal_p", prefix="metal_p_test"
            )
            assert cmd.count_atoms(atoms["donor_p"]) == 1
            assert cmd.count_atoms(atoms["donors"]) == 1
        finally:
            cmd.delete("all")

    def test_frame_swallows_origin_and_rebuild_exceptions(self, mocker):
        """``cmd.origin``/``cmd.rebuild`` failures are swallowed so the
        rest of ``frame()`` still runs (older PyMOL builds lack these)."""
        mock_cmd = mocker.patch(
            "chemsmart.jobs.mol.templates.zhang_group_scientific_styles.cmd"
        )
        mock_cmd.origin.side_effect = Exception("no origin support")
        mock_cmd.rebuild.side_effect = Exception("no rebuild support")

        style = ScientificStyle()
        style.frame("all", "core_name")

        mock_cmd.origin.assert_called_once_with("core_name")
        mock_cmd.rebuild.assert_called_once()
        mock_cmd.refresh.assert_called_once()


@pytest.mark.usefixtures("skip_if_no_pymol")
class TestPyMOLScientificStyleDefensiveExceptionBranches:
    """Force individual real-PyMOL API calls to fail (via targeted
    ``mocker.patch.object`` on the actual ``pymol.cmd`` module) to exercise
    this file's ``except Exception: pass`` guards, which real PyMOL calls
    practically never trigger under normal, valid-selection usage."""

    def test_safe_ray_shadows_swallows_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(
            cmd.util, "ray_shadows", side_effect=Exception("boom")
        )
        ScientificStyle.safe_ray_shadows("light")

    def test_hide_distance_value_labels_swallows_hide_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(
            ScientificStyle, "_distance_object_names", return_value=["d1"]
        )
        mocker.patch.object(cmd, "hide", side_effect=Exception("boom"))
        ScientificStyle.hide_distance_value_labels()

    def test_distance_object_names_swallows_get_names_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(cmd, "get_names", side_effect=Exception("boom"))
        assert ScientificStyle._distance_object_names() == []

    def test_pairs_from_distance_objects_swallows_get_state_exception(
        self, mocker
    ):
        from pymol import cmd

        mocker.patch.object(cmd, "get_state", side_effect=Exception("boom"))
        assert ScientificStyle.pairs_from_distance_objects() == []

    def test_pairs_from_distance_objects_swallows_iterate_state_exception(
        self, mocker
    ):
        from pymol import cmd

        mocker.patch.object(
            ScientificStyle, "_distance_object_names", return_value=["d1"]
        )
        mocker.patch.object(
            cmd, "iterate_state", side_effect=Exception("boom")
        )
        assert ScientificStyle.pairs_from_distance_objects() == []

    def test_pairs_from_distance_objects_swallows_get_session_exception(
        self, mocker
    ):
        from pymol import cmd

        mocker.patch.object(
            ScientificStyle, "_distance_object_names", return_value=["d1"]
        )
        mocker.patch.object(cmd, "get_session", side_effect=Exception("boom"))
        assert ScientificStyle.pairs_from_distance_objects() == []

    def test_bond_atom_index_pairs_swallows_bond_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(cmd, "bond", side_effect=Exception("boom"))
        ScientificStyle.bond_atom_index_pairs([(1, 2)])

    def test_remove_distance_objects_swallows_delete_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(
            ScientificStyle, "_distance_object_names", return_value=["d1"]
        )
        mocker.patch.object(cmd, "delete", side_effect=Exception("boom"))
        ScientificStyle.remove_distance_objects()

    def test_apply_element_palette_swallows_color_exceptions(self, mocker):
        from pymol import cmd

        mocker.patch.object(cmd, "color", side_effect=Exception("boom"))
        ScientificStyle.apply_element_palette(
            "all", {"C": "black"}, overrides={"some_sel": "gold"}
        )

    def test_define_colors_swallows_set_color_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(cmd, "set_color", side_effect=Exception("boom"))
        ComicMetallicStyle().define_colors()

    def test_define_shared_colors_swallows_set_color_exception(self, mocker):
        from pymol import cmd

        mocker.patch.object(cmd, "set_color", side_effect=Exception("boom"))
        ScientificStyle.define_shared_colors()

    def test_soft_cartoon_render_swallows_set_bond_exception(
        self, loaded_1_mer, mocker
    ):
        from pymol import cmd

        mocker.patch.object(cmd, "set_bond", side_effect=Exception("boom"))
        SoftCartoonStyle().render(loaded_1_mer)

    def test_neon_coordination_core_render_swallows_set_bond_exception(
        self, loaded_1_mer, mocker
    ):
        from pymol import cmd

        mocker.patch.object(cmd, "set_bond", side_effect=Exception("boom"))
        NeonCoordinationCoreStyle().render(loaded_1_mer)

    def test_quasi_chemdraw_bold_render_swallows_set_bond_exception(
        self, loaded_1_mer, mocker
    ):
        from pymol import cmd

        mocker.patch.object(cmd, "set_bond", side_effect=Exception("boom"))
        QuasiChemDrawBoldStyle().render(loaded_1_mer)

    def test_pairs_from_distance_objects_skips_none_and_malformed_and_unmatched_points(
        self, mocker
    ):
        """A ``d1``/``d2``/… session's per-object point payload can be
        ``None`` (line-object with no visible dash points), malformed
        (older/newer PyMOL session dict layouts), or reference coordinates
        no longer present in ``xyz2idx`` (selection changed since the
        distance was drawn). All three must be skipped, not raised."""
        from pymol import cmd

        mocker.patch.object(
            ScientificStyle,
            "_distance_object_names",
            return_value=["d1", "d2", "d3"],
        )
        mocker.patch.object(cmd, "get_state", return_value=1)
        mocker.patch.object(cmd, "iterate_state", return_value=None)

        def make_session_object(points):
            # obj[5][2][state - 1][1] == points, with state == 1.
            inner = [None, points]
            return [None] * 5 + [[None, None, [inner]]]

        none_points_obj = make_session_object(None)
        malformed_obj = [None] * 3  # obj[5] raises IndexError
        unmatched_obj = make_session_object(
            [10.0, 10.0, 10.0, 20.0, 20.0, 20.0]
        )
        mocker.patch.object(
            cmd,
            "get_session",
            return_value={
                "names": [none_points_obj, malformed_obj, unmatched_obj]
            },
        )

        # xyz2idx stays empty since `iterate_state` is mocked out, so the
        # unmatched_obj's points can never resolve to an atom index.
        assert ScientificStyle.pairs_from_distance_objects("all") == []
