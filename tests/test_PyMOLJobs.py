import inspect
import os.path
import shutil
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
    PyMOLJobRunner,
    PyMOLNCIJobRunner,
    PyMOLScientificStyleVisualizationJobRunner,
    PyMOLSpinJobRunner,
    normalize_pymol_style,
)
from chemsmart.jobs.mol.spin import PyMOLSpinJob
from chemsmart.jobs.mol.templates.scientific_styles import (
    ELEMENT_CATEGORIES,
    ComicMetallicStyle,
    EditorialMinimalStyle,
    MatteClayStyle,
    NeonCoordinationCoreStyle,
    QuasiChemDrawBoldStyle,
    ScientificStyle,
    SoftCartoonStyle,
    SoftCeramicStyle,
    StericSurfaceStyle,
    XrayWireStyle,
    _get_coordinating_atoms,
    element_category_selection,
    hide_distance_value_labels,
    render_editorial_minimal,
    render_matte_clay,
    render_neon_coordination_core,
    render_quasi_chemdraw_bold,
    render_soft_cartoon,
    render_soft_ceramic,
    render_steric_surface,
    render_xray_wire,
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

    def test_1_mer_xyz_is_valid_visualization_geometry(
        self, visualized_1_mer_xyz_file
    ):
        """``1-mer.xyz`` loads as a single Mn coordination complex."""
        molecules = Molecule.from_filepath(
            visualized_1_mer_xyz_file, return_list=True
        )

        assert len(molecules) == 1
        assert molecules[0].num_atoms == 36
        assert "Mn" in molecules[0].elements

    def test_format_pymol_style_command_is_independent_of_coordinates(self):
        """``-c`` is handled via distance/angle labels, not style args."""
        for style, expected in (
            (
                "comic",
                f"render_comic_metallic_labeled_final {self.label_1_mer}",
            ),
            ("soft_cartoon", f"render_soft_cartoon {self.label_1_mer}"),
            (
                "neon_coordination_core",
                f"render_neon_coordination_core {self.label_1_mer}",
            ),
            (
                "editorial_minimal",
                f"render_editorial_minimal {self.label_1_mer}",
            ),
            ("soft_ceramic", f"render_soft_ceramic {self.label_1_mer}"),
            ("matte_clay", f"render_matte_clay {self.label_1_mer}"),
            (
                "glossy",
                f"metallic_poster_render {self.label_1_mer}, elem Mn, None, "
                f"2.6, N+O+S+P+H",
            ),
        ):
            for coordinates in (None, self.coordination_bonds_1_mer):
                job = SimpleNamespace(style=style, coordinates=coordinates)
                command = PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                    job, self.label_1_mer
                )
                assert command == expected

    def test_scientific_styles_hide_distance_value_labels(self):
        """Derived styles keep distance dashes but hide numeric distance labels."""
        hide_source = inspect.getsource(hide_distance_value_labels)
        assert 'cmd.hide("labels"' in hide_source
        assert "suffix.isdigit()" in hide_source

        finish_source = inspect.getsource(ScientificStyle.finish_default)
        assert "self.finalize()" in finish_source

        comic_source = inspect.getsource(ComicMetallicStyle.render)
        assert "self.finalize()" in comic_source

    def test_scientific_style_runner_orders_distances_before_render(self):
        """``-c`` distances must exist before style render calls finalize."""
        source = inspect.getsource(
            PyMOLScientificStyleVisualizationJobRunner._get_command
        )
        assert source.index("_add_coordinates_labels") < source.index(
            "_append_style_render"
        )
        setup_source = inspect.getsource(
            PyMOLScientificStyleVisualizationJobRunner._setup_style
        )
        assert ' -d "' in setup_source
        assert "render_" not in setup_source

    def test_get_coordinating_atoms_uses_radius_ratio_helper(self):
        source = inspect.getsource(_get_coordinating_atoms)
        module_source = inspect.getsource(
            inspect.getmodule(_get_coordinating_atoms)
        )

        assert (
            "from chemsmart.utils.geometry import get_coordinating_atoms"
            in (module_source)
        )
        assert "get_coordinating_atoms(" in source
        assert "cmd.get_model" in source
        assert "primary_local" in source
        assert "secondary_local" in source
        assert "coordination_core" in source
        assert "within %s of" not in source
        assert "cmd.color" not in source
        assert "sphere_scale" not in source

    def test_element_category_helpers_are_centralized(self):
        assert ELEMENT_CATEGORIES["halogen"] == "F+Cl+Br+I"
        assert ELEMENT_CATEGORIES["N+O"] == "N+O"
        assert ELEMENT_CATEGORIES["metal"].startswith("elem ")
        assert (
            element_category_selection("all", "halogen")
            == "(all) and elem F+Cl+Br+I"
        )
        assert (
            element_category_selection("sc_shell", "N+O")
            == "(sc_shell) and elem N+O"
        )

        base_source = inspect.getsource(ScientificStyle)
        assert (
            "def _safe_set(self, setting, value, selection=None, category=None)"
            in base_source
        )
        assert (
            "def apply_style_palette(self, selection, palette, overrides=None)"
            in base_source
        )
        assert "element_category_selection" in base_source

    def test_render_editorial_minimal_defines_expected_visual_parameters(self):
        source = inspect.getsource(EditorialMinimalStyle)
        wrapper = inspect.getsource(render_editorial_minimal)

        assert issubclass(EditorialMinimalStyle, ScientificStyle)
        assert EditorialMinimalStyle.prefix == "editorial"
        assert "EditorialMinimalStyle().render" in wrapper
        assert "select_coordination" in source
        assert 'metal="elem Mn"' not in source
        assert "donors=" not in source
        assert "within 2.8" not in source
        assert '_safe_set("sphere_scale", 0.60, atoms["metal"])' in source
        assert '_safe_set("sphere_scale", 0.36, atoms["donor_n"])' in source
        assert '_safe_set("sphere_scale", 0.39, atoms["donor_s"])' in source
        assert '_safe_set("sphere_scale", 0.26, atoms["co_c"])' in source
        assert '_safe_set("sphere_scale", 0.21, atoms["co_o"])' in source
        assert '_safe_set("stick_radius", 0.12, sel)' in source
        assert (
            '_safe_set("stick_radius", 0.15, atoms["coordination_core"])'
            in source
        )
        assert '_safe_set("field_of_view", 45)' in source
        assert 'cmd.color("mn_rose", atoms["metal"])' in source
        assert 'cmd.color("sulfur_gold", atoms["donor_s"])' in source
        assert '_safe_set("ambient_occlusion_mode", 1)' in source

    def test_render_neon_coordination_core_defines_expected_visual_parameters(
        self,
    ):
        source = inspect.getsource(NeonCoordinationCoreStyle)
        wrapper = inspect.getsource(render_neon_coordination_core)

        assert issubclass(NeonCoordinationCoreStyle, ScientificStyle)
        assert NeonCoordinationCoreStyle.prefix == "ncc"
        assert "NeonCoordinationCoreStyle().render" in wrapper
        assert "select_coordination" in source
        assert "_common_select_core" not in source
        assert "within 3.00" not in source
        assert "neighbor ncc_metal" not in source
        assert '_safe_set("sphere_scale", 0.44, atoms["metal"])' in source
        assert '_safe_set("stick_radius", 0.10, sel)' in source
        assert 'cmd.set_bond("stick_radius", 0.145' in source
        assert "apply_style_palette" in source
        assert 'overrides={atoms["metal"]: "ncc_metal_c"}' in source
        assert '_safe_set("opaque_background", 0)' in source
        assert '_safe_set("ray_opaque_background", 0)' in source
        assert '_safe_set("ambient_occlusion_mode", 1)' in source
        assert '_safe_set("depth_cue", 0)' in source

    def test_render_soft_ceramic_defines_expected_visual_parameters(self):
        source = inspect.getsource(SoftCeramicStyle)
        wrapper = inspect.getsource(render_soft_ceramic)

        assert issubclass(SoftCeramicStyle, ScientificStyle)
        assert SoftCeramicStyle.prefix == "soft_ceramic"
        assert SoftCeramicStyle.include_nh_h is True
        assert "SoftCeramicStyle().render" in wrapper
        assert "select_coordination" in source
        assert 'metal="elem Mn"' not in source
        assert "donors=" not in source
        assert '_safe_set("sphere_scale", 0.56, atoms["metal"])' in source
        assert '_safe_set("sphere_scale", 0.40, atoms["donor_s"])' in source
        assert '_safe_set("sphere_scale", 0.37, atoms["donor_n"])' in source
        assert '_safe_set("sphere_scale", 0.28, atoms["co_c"])' in source
        assert '_safe_set("sphere_scale", 0.18, atoms["co_o"])' in source
        assert '_safe_set("sphere_scale", 0.25, atoms["hydride"])' in source
        assert '_safe_set("stick_radius", 0.115, sel)' in source
        assert 'cmd.color("metal_bronze", atoms["metal"])' in source
        assert 'cmd.color("sulfur_soft_gold", atoms["donor_s"])' in source
        assert 'cmd.bg_color("studio_background")' in source
        assert '_safe_set("ray_opaque_background", 1)' in source

    def test_render_soft_cartoon_defines_expected_visual_parameters(self):
        source = inspect.getsource(SoftCartoonStyle)
        wrapper = inspect.getsource(render_soft_cartoon)

        assert issubclass(SoftCartoonStyle, ScientificStyle)
        assert SoftCartoonStyle.prefix == "sc"
        assert SoftCartoonStyle.include_nh_h is True
        assert "SoftCartoonStyle().render" in wrapper
        assert "select_coordination" in source
        assert "within 2.6" not in source
        assert "byres" not in source
        assert "elem Mn or elem Fe" not in source
        assert 'cmd.bg_color("sc_background")' in source
        assert '_safe_set("stick_radius", 0.135, sel)' in source
        assert 'cmd.set_bond("stick_radius", 0.165' in source
        assert '_safe_set("sphere_scale", 0.40, atoms["metal"])' in source
        assert (
            '_safe_set("sphere_scale", 0.215, f"({shell}) and not elem H")'
            in source
        )
        assert '("N+O", 0.275)' in source
        assert "category=category" in source
        assert "self._safe_set" in source
        assert "apply_style_palette" in source
        assert 'overrides={atoms["metal"]: "sc_metal"}' in source
        assert "_apply_lighting(" in source
        assert "apply_soft_shadows" in source
        assert "apply_ambient_occlusion" in source
        assert "apply_illustrated_camera" in source
        assert '_safe_set("ray_trace_mode", 1)' in source
        assert "self.frame(" in source

    def test_render_quasi_chemdraw_bold_defines_expected_visual_parameters(
        self,
    ):
        style_source = inspect.getsource(QuasiChemDrawBoldStyle)
        wrapper_source = inspect.getsource(render_quasi_chemdraw_bold)

        assert issubclass(QuasiChemDrawBoldStyle, ScientificStyle)
        assert QuasiChemDrawBoldStyle.prefix == "qcd"
        assert QuasiChemDrawBoldStyle.command == "render_quasi_chemdraw_bold"
        assert "qcd_bond" in QuasiChemDrawBoldStyle.colors
        assert "QuasiChemDrawBoldStyle().render" in wrapper_source
        assert "select_coordination" in style_source
        assert "_common_select_core" not in style_source
        assert "within 2.75" not in style_source
        assert "neighbor qcd_metal" not in style_source
        assert "Li+Na+K+Rb+Cs" not in style_source
        assert '_safe_set("stick_radius", 0.155, sel)' in style_source
        assert '_safe_set("stick_color", "qcd_bond", sel)' in style_source
        assert 'cmd.set_bond("stick_radius", 0.185' in style_source
        assert (
            '_safe_set("sphere_scale", 0.34, atoms["metal"])' in style_source
        )
        assert "apply_style_palette" in style_source
        assert 'overrides={atoms["metal"]: "qcd_metal"}' in style_source
        assert '_safe_set("ambient", 0.58)' in style_source
        assert '_safe_set("ray_shadow", 0)' in style_source
        assert '_safe_set("ambient_occlusion_mode", 0)' in style_source
        assert '_safe_set("orthoscopic", 1)' in style_source

    def test_render_xray_wire_defines_expected_visual_parameters(self):
        style_source = inspect.getsource(XrayWireStyle)
        wrapper_source = inspect.getsource(render_xray_wire)

        assert issubclass(XrayWireStyle, ScientificStyle)
        assert XrayWireStyle.prefix == "xray"
        assert XrayWireStyle.command == "render_xray_wire"
        assert "xw_charcoal" in XrayWireStyle.colors
        assert "xw_metal" in XrayWireStyle.colors
        assert "XrayWireStyle().render" in wrapper_source
        assert "select_coordination" in style_source
        assert "_begin_scientific_style" not in style_source
        assert "_common_select_core" not in style_source
        assert "_finish_scientific_style" not in style_source
        assert 'cmd.label("all"' not in style_source
        assert "id 0" not in style_source
        assert 'center_metal = atoms["metal"]' in style_source
        assert 'cmd.hide("lines", sel)' in style_source
        assert 'cmd.hide("nonbonded", sel)' in style_source
        assert '_safe_set("stick_radius", 0.18, sel)' in style_source
        assert '_safe_set("stick_ball", 1)' in style_source
        assert '_safe_set("stick_ball_ratio", 1.0)' in style_source
        assert '_safe_set("sphere_scale", 0.6, center_metal)' in style_source
        assert "apply_style_palette" in style_source
        assert 'overrides={center_metal: "xw_metal"}' in style_source
        assert '_safe_set("ray_trace_mode", 3)' in style_source
        assert '_safe_set("ray_trace_gain", 0.05)' in style_source
        assert '_safe_set("ambient", 0.5)' in style_source
        assert '_safe_set("direct", 0.6)' in style_source
        assert '_safe_set("reflect", 0.0)' in style_source
        assert '_safe_set("spec_power", 1.0)' in style_source
        assert '_safe_set("spec_count", 0)' in style_source
        assert '_safe_set("ray_shadow", 1)' in style_source
        assert '_safe_set("ray_trace_fog", 0)' in style_source
        assert "_set_transparent_background()" in style_source
        assert "apply_illustrated_camera" in style_source
        assert "self.finish_default(selection)" in style_source

    def test_render_steric_surface_defines_expected_visual_parameters(self):
        style_source = inspect.getsource(StericSurfaceStyle)
        wrapper_source = inspect.getsource(render_steric_surface)

        assert issubclass(StericSurfaceStyle, ScientificStyle)
        assert StericSurfaceStyle.prefix == "steric"
        assert StericSurfaceStyle.command == "render_steric_surface"
        assert "StericSurfaceStyle().render" in wrapper_source
        assert "select_coordination" in style_source
        assert "_begin_scientific_style" not in style_source
        assert "_common_select_core" not in style_source
        assert "_show_coord_ball_and_stick" not in style_source
        assert "_apply_coord_sphere_scales" not in style_source
        assert "_color_by_element" not in style_source
        assert "_finish_scientific_style" not in style_source
        assert "within 2.6" not in style_source
        assert "byres" not in style_source
        assert "sphere_atoms" in style_source
        assert 'cmd.show("spheres", sphere_atoms)' in style_source
        assert (
            'cmd.show("spheres", atoms["coordination_core"])'
            not in style_source
        )
        assert (
            '_safe_set("sphere_scale", 0.52, atoms["metal"])' in style_source
        )
        assert '_safe_set("sphere_scale", 0.34, atoms["co_c"])' in style_source
        assert '_safe_set("sphere_scale", 0.34, atoms["co_o"])' in style_source
        assert 'cmd.color("sci_C_gray", atoms["co_c"])' in style_source
        assert 'cmd.color("sci_O_red", atoms["co_o"])' in style_source
        assert "apply_style_palette" in style_source
        assert 'overrides={atoms["metal"]: "metal_gold"}' in style_source
        assert '_safe_set("transparency", 0.68, sel)' in style_source
        assert "self.finish_default(selection)" in style_source

    def test_render_matte_clay_defines_expected_visual_parameters(self):
        source = inspect.getsource(MatteClayStyle)
        wrapper = inspect.getsource(render_matte_clay)

        assert issubclass(MatteClayStyle, ScientificStyle)
        assert MatteClayStyle.prefix == "matte_clay"
        assert "MatteClayStyle().render" in wrapper
        assert "select_coordination" in source
        assert "donors=" not in source
        assert 'atoms["hydride"]' in source
        assert '_safe_set("sphere_scale", 0.62, atoms["metal"])' in source
        assert 'cmd.color("mc_metal_center", atoms["metal"])' in source
        assert 'cmd.bg_color("white")' in source
        assert '_safe_set("ray_opaque_background", 0)' in source
        assert '_safe_set("ambient_occlusion_mode", 1)' in source
        assert '_safe_set("ambient_occlusion_scale", 18)' in source
        assert "cmd.select(" not in source

    def test_hybrid_is_not_a_derived_style(self):
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
        assert "hide_distance_value_labels" in inspect.getsource(
            hide_distance_value_labels
        )


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
        assert os.path.exists(os.path.join(tmpdir, "scientific_styles.py"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"render_comic_metallic_labeled_final {job.label}"
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
            == f"render_comic_metallic_labeled_final {job.label}"
        )
        runner = PyMOLScientificStyleVisualizationJobRunner.__new__(
            PyMOLScientificStyleVisualizationJobRunner
        )
        assert "distance d1, id 1, id 2" in runner._add_coordinates_labels(
            job, "cmd"
        )
        assert "hide labels" not in runner._add_coordinates_labels(job, "cmd")
        assert "self.finalize()" in inspect.getsource(
            ComicMetallicStyle.render
        )

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
        assert os.path.exists(os.path.join(tmpdir, "scientific_styles.py"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            ).startswith("metallic_poster_render")
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
        assert os.path.exists(os.path.join(tmpdir, "scientific_styles.py"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"render_editorial_minimal {job.label}"
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
        assert os.path.exists(os.path.join(tmpdir, "scientific_styles.py"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.xyz"))
        assert os.path.exists(os.path.join(tmpdir, f"{job.label}.pse"))
        assert (
            PyMOLScientificStyleVisualizationJobRunner._format_style_command(
                job, job.label
            )
            == f"render_soft_ceramic {job.label}"
        )
