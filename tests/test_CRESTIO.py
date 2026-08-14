import os

import numpy as np

from chemsmart.io.crest.file import CRESTEnergiesFile, CRESTMainOut
from chemsmart.io.crest.folder import CRESTFolder
from chemsmart.io.crest.output import CRESTOutput
from chemsmart.io.crest.route import CRESTRoute
from chemsmart.utils.constants import kcal_per_mol_to_hartree


class TestCRESTRoute:
    """Tests for CRESTRoute class."""

    def test_read_route(self):
        s0 = "crest mol.xyz"
        r0 = CRESTRoute(route_string=s0)
        assert r0.gfn_version is None
        assert r0.optimization_level is None
        assert r0.solvent_model is None
        assert r0.solvent_id is None
        assert r0.energy_window is None
        assert r0.quick_mode is None
        assert r0.threads is None
        assert r0.nci is False
        assert r0.constrained is False

        s1 = "crest mol.xyz --gfn2//gfnff"
        r1 = CRESTRoute(route_string=s1)
        assert r1.method == "gfn2//gfnff"
        assert r1.gfn_version == "gfn2//gfnff"
        assert r1.basis == "default"
        assert r1.charge == 0
        assert r1.uhf == 0
        assert r1.multiplicity == 1
        assert r1.jobtype == "conformers"

        s2 = (
            "crest mol.xyz "
            "--gfn2 "
            "--chrg -1 "
            "--uhf 2 "
            "--alpb water "
            "--optlev tight "
            "--ewin 30 "
            "--squick "
            "-T 64"
        )
        r2 = CRESTRoute(route_string=s2)
        assert r2.method == "gfn2"
        assert r2.gfn_version == "gfn2"
        assert r2.charge == -1
        assert r2.uhf == 2
        assert r2.multiplicity == 3
        assert r2.solvent_model == "alpb"
        assert r2.solvent_id == "water"
        assert r2.optimization_level == "tight"
        assert r2.energy_window == 30.0
        assert r2.quick_mode == "squick"
        assert r2.threads == 64

        # GBSA solvent
        s3 = "crest mol.xyz --gfn2 --gbsa toluene"
        r3 = CRESTRoute(route_string=s3)
        assert r3.solvent_model == "gbsa"
        assert r3.solvent_id == "toluene"

        # CREGEN parameters
        s4 = (
            "crest mol.xyz "
            "--gfn2 "
            "--rthr 0.15 "
            "--ethr 0.10 "
            "--bthr 0.02 "
            "--pthr 0.10"
        )
        r4 = CRESTRoute(route_string=s4)
        assert r4.rmsd_threshold == 0.15
        assert r4.energy_threshold == 0.10
        assert r4.bconst_threshold == 0.02
        assert r4.population_threshold == 0.10

        # MD options
        s5 = (
            "crest mol.xyz "
            "--gfn2 "
            "--temp 300 "
            "--shake 1 "
            "--tstep 2 "
            "--mdlen 50 "
            "--mddump 200 "
            "--vbdump 0.5 "
            "--tnmd 500"
        )
        r5 = CRESTRoute(route_string=s5)
        assert r5.temperature == 300.0
        assert r5.shake == 1
        assert r5.md_timestep == 2
        assert r5.md_length == 50.0
        assert r5.md_dump_step == 200
        assert r5.vbias_dump_interval == 0.5
        assert r5.additional_md_temperature == 500.0

        # Scaled MD length
        s6 = "crest mol.xyz --gfn2 -len x2.0"
        r6 = CRESTRoute(route_string=s6)
        assert r6.md_length == "x2.0"

        # Boolean flags
        s7 = (
            "crest mol.xyz "
            "--gfn2 "
            "--nci "
            "--cinp constraints.inp "
            "--noreftopo "
            "--notopo"
        )
        r7 = CRESTRoute(route_string=s7)
        assert r7.nci is True
        assert r7.constrained is True
        assert r7.no_reference_topology_check is True
        assert r7.no_topology_check is True

        s8 = "crest mol.xyz --gff --opt verytight"
        r8 = CRESTRoute(route_string=s8)
        assert r8.gfn_version == "gfnff"
        assert r8.optimization_level == "vtight"

        s9 = (
            "crest mol.xyz --gfn2 "
            "-nci -cinp constraints.inp -noreftopo -notopo "
            "-ewin 10 -rthr 0.2 -ethr 0.05 -bthr 0.01 -pthr 0.1 "
            "-temp 250 -shake 0 -tstep 3 -mdlen 20 "
            "-mddump 50 -vbdump 2.0 -tnmd 350"
        )
        r9 = CRESTRoute(route_string=s9)
        assert r9.nci is True
        assert r9.constrained is True
        assert r9.no_reference_topology_check is True
        assert r9.no_topology_check is True
        assert r9.energy_window == 10.0
        assert r9.rmsd_threshold == 0.2
        assert r9.energy_threshold == 0.05
        assert r9.bconst_threshold == 0.01
        assert r9.population_threshold == 0.1
        assert r9.temperature == 250.0
        assert r9.shake == 0
        assert r9.md_timestep == 3
        assert r9.md_length == 20.0
        assert r9.md_dump_step == 50
        assert r9.vbias_dump_interval == 2.0
        assert r9.additional_md_temperature == 350.0


class TestCRESTMainOut:
    """Tests for XTBMainOut class."""

    def test_main_out_octane(self, crest_octane_outfolder):
        """Test parsing main output from octane conformational sampling."""
        octane_main_out_file = os.path.join(
            crest_octane_outfolder, "octane_conformers.out"
        )
        assert os.path.exists(octane_main_out_file)
        octane_main_out = CRESTMainOut(octane_main_out_file)
        assert octane_main_out.version == "3.0.2"
        assert octane_main_out.normal_termination
        assert (
            octane_main_out.route_string
            == "crest octane.xyz --gfn2 --chrg 0 --uhf 0 --optlev vtight -T 64"
        )
        assert not octane_main_out.constrained
        assert not octane_main_out.topology_mismatch
        assert octane_main_out.num_conformers == 104
        assert not octane_main_out.single_conformer
        assert octane_main_out.num_rotamers == 1342
        assert octane_main_out.num_atoms == 26
        assert octane_main_out.rmsd_threshold == 0.125
        assert octane_main_out.bconst_threshold == 0.01
        assert octane_main_out.population_threshold == 0.05
        assert octane_main_out.energy_window == 6.0
        assert np.isclose(octane_main_out.reference_state_energy, -26.32263439)
        assert octane_main_out.temperature == 298.15
        assert np.isclose(octane_main_out.lowest_energy, -26.32263)
        relative_average_energy_ensemble_kcal_per_mol = (
            1.001  # <E_rel> in kcal/mol
        )
        entropy_ensemble_cal_per_mol_K = 12.422  # S_ensemble in cal/mol K
        minus_TS_kcal_per_mol = -3.703  # -T*S_ensemble in kcal/mol
        # E_ensemble = E_reference + <E_rel>
        expected_average_energy_ensemble = (
            octane_main_out.reference_state_energy
            + relative_average_energy_ensemble_kcal_per_mol
            * kcal_per_mol_to_hartree
        )
        assert np.isclose(
            octane_main_out.average_energy_ensemble,
            expected_average_energy_ensemble,
        )
        expected_entropy_ensemble = (
            entropy_ensemble_cal_per_mol_K * 1e-3 * kcal_per_mol_to_hartree
        )
        assert np.isclose(
            octane_main_out.entropy_ensemble, expected_entropy_ensemble
        )
        # G_ensemble = E_reference + <E_rel> - T*S_ensemble
        expected_free_energy_ensemble = (
            octane_main_out.reference_state_energy
            + (
                relative_average_energy_ensemble_kcal_per_mol
                + minus_TS_kcal_per_mol
            )
            * kcal_per_mol_to_hartree
        )
        assert np.isclose(
            octane_main_out.free_energy_ensemble, expected_free_energy_ensemble
        )
        assert np.isclose(
            octane_main_out.free_energy_ensemble,
            octane_main_out.average_energy_ensemble
            - octane_main_out.temperature * octane_main_out.entropy_ensemble,
        )
        assert np.isclose(octane_main_out.population_of_lowest, 8.348)
        # ~36 min 50 s ≈ 0.614 hours
        assert np.isclose(octane_main_out.wall_time, 0.614, atol=1e-2)
        # ~14 h 48 min 33 s ≈ 14.809
        assert np.isclose(octane_main_out.cpu_time, 14.809, atol=1e-2)

    def test_main_out_ts1a_constrained(self, crest_ts1a_constrained_outfolder):
        """Test parsing main output from TS1A constrained sampling."""
        ts1a_main_out_file = os.path.join(
            crest_ts1a_constrained_outfolder, "TS1A_conformers.out"
        )
        assert os.path.exists(ts1a_main_out_file)
        ts1a_main_out = CRESTMainOut(ts1a_main_out_file)
        assert ts1a_main_out.version == "3.0.2"
        assert ts1a_main_out.normal_termination
        assert (
            ts1a_main_out.route_string
            == "crest TS1A.xyz --cinp constraints.inp --gfn2 --chrg 0 --uhf 0 --optlev tight -T 64"
        )
        assert ts1a_main_out.constrained
        assert ts1a_main_out.num_conformers == 22
        assert ts1a_main_out.num_rotamers == 63
        assert ts1a_main_out.num_atoms == 40

    def test_main_out_hexane_failed(self, crest_hexane_failed_outfolder):
        """Failed CREST run: parseable header/route, but no ensemble summary."""
        hexane_main_out_file = os.path.join(
            crest_hexane_failed_outfolder, "hexane.out"
        )
        assert os.path.exists(hexane_main_out_file)
        hexane_main_out = CRESTMainOut(hexane_main_out_file)
        assert hexane_main_out.version == "3.0.2"
        assert not hexane_main_out.normal_termination
        assert hexane_main_out.route_string.startswith("crest ")
        assert "hexane.xyz" in hexane_main_out.route_string
        assert hexane_main_out.route_string.endswith(
            "--gfn2 --chrg 0 --uhf 0 --ewin 6.0 --optlev vtight"
        )
        assert not hexane_main_out.constrained
        assert not hexane_main_out.topology_mismatch
        assert hexane_main_out.charge == 0
        assert hexane_main_out.uhf == 0
        assert hexane_main_out.multiplicity == 1
        assert hexane_main_out.num_conformers is None
        assert hexane_main_out.num_rotamers is None
        assert hexane_main_out.num_atoms is None
        assert hexane_main_out.lowest_energy is None
        assert hexane_main_out.wall_time is None
        assert hexane_main_out.cpu_time is None


class TestCRESTEnergiesFile:
    """Tests for CRESTEnergiesFile class."""

    def test_energies_octane(self, crest_octane_outfolder):
        """Test parsing octane energies from CREST energies file."""
        octane_energies_file = os.path.join(
            crest_octane_outfolder, "crest.energies"
        )
        octane_energies = CRESTEnergiesFile(octane_energies_file)
        assert octane_energies.num_conformers == 104
        assert np.isclose(octane_energies.relative_energies[0], 0.000)
        assert np.isclose(octane_energies.relative_energies[1], 0.544)
        assert np.isclose(octane_energies.relative_energies[-1], 5.404)
        assert np.isclose(octane_energies.energy_range, 5.404 - 0.000)

    def test_energies_ts1a_constrained(self, crest_ts1a_constrained_outfolder):
        """Test parsing ts1a energies from CREST energies file."""
        ts1a_energies_file = os.path.join(
            crest_ts1a_constrained_outfolder, "crest.energies"
        )
        ts1a_energies = CRESTEnergiesFile(ts1a_energies_file)
        assert ts1a_energies.num_conformers == 22
        assert np.isclose(ts1a_energies.relative_energies[0], 0.000)
        assert np.isclose(ts1a_energies.relative_energies[1], 2.901)
        assert np.isclose(ts1a_energies.relative_energies[-2], 4.041)
        assert np.isclose(ts1a_energies.relative_energies[-1], 4.031)
        assert np.isclose(ts1a_energies.energy_range, 4.041 - 0.000)

    def test_energies_styrene_single(self, crest_styrene_outfolder):
        """Test parsing styrene energies from CREST energies file."""
        styrene_energies_file = os.path.join(
            crest_styrene_outfolder, "crest.energies"
        )
        styrene_energies = CRESTEnergiesFile(styrene_energies_file)
        assert styrene_energies.num_conformers == 1
        assert np.isclose(styrene_energies.relative_energies[0], 0.000)
        assert np.isclose(styrene_energies.energy_range, 0.000)


class TestCRESTFolder:
    """Tests for CRESTFolder class."""

    def test_folder_octane(self, crest_octane_outfolder):
        """Test CRESTFolder class with octane conformational sampling output."""
        assert os.path.exists(crest_octane_outfolder)
        octane_folder = CRESTFolder(crest_octane_outfolder)
        assert octane_folder.is_crest_calculation_directory

        assert octane_folder.crest_out_filepath is not None
        assert (
            os.path.basename(octane_folder.crest_out_filepath)
            == "octane_conformers.out"
        )

        assert octane_folder.conformers_xyz_filepath is not None
        assert (
            os.path.basename(octane_folder.conformers_xyz_filepath)
            == "crest_conformers.xyz"
        )

        assert octane_folder.best_xyz_filepath is not None
        assert (
            os.path.basename(octane_folder.best_xyz_filepath)
            == "crest_best.xyz"
        )

        assert octane_folder.rotamers_xyz_filepath is not None
        assert (
            os.path.basename(octane_folder.rotamers_xyz_filepath)
            == "crest_rotamers.xyz"
        )

        assert octane_folder.energies_filepath is not None
        assert (
            os.path.basename(octane_folder.energies_filepath)
            == "crest.energies"
        )

        assert octane_folder.constraints_inp_filepath is None

        assert octane_folder.crestopt_log_filepath is not None
        assert (
            os.path.basename(octane_folder.crestopt_log_filepath)
            == "crestopt.log"
        )

        assert octane_folder.coord_filepath is not None
        assert os.path.basename(octane_folder.coord_filepath) == "coord"

    def test_folder_ts1a_constrained(self, crest_ts1a_constrained_outfolder):
        """Test CRESTFolder class with TS1A constrained sampling output."""
        assert os.path.exists(crest_ts1a_constrained_outfolder)
        ts1a_folder = CRESTFolder(crest_ts1a_constrained_outfolder)
        assert ts1a_folder.is_crest_calculation_directory

        assert ts1a_folder.crest_out_filepath is not None
        assert (
            os.path.basename(ts1a_folder.crest_out_filepath)
            == "TS1A_conformers.out"
        )

        assert ts1a_folder.conformers_xyz_filepath is not None
        assert (
            os.path.basename(ts1a_folder.conformers_xyz_filepath)
            == "crest_conformers.xyz"
        )

        assert ts1a_folder.best_xyz_filepath is not None
        assert (
            os.path.basename(ts1a_folder.best_xyz_filepath) == "crest_best.xyz"
        )

        assert ts1a_folder.rotamers_xyz_filepath is not None
        assert (
            os.path.basename(ts1a_folder.rotamers_xyz_filepath)
            == "crest_rotamers.xyz"
        )

        assert ts1a_folder.energies_filepath is not None
        assert (
            os.path.basename(ts1a_folder.energies_filepath) == "crest.energies"
        )

        assert ts1a_folder.constraints_inp_filepath is not None
        assert (
            os.path.basename(ts1a_folder.constraints_inp_filepath)
            == "constraints.inp"
        )

        assert ts1a_folder.crestopt_log_filepath is not None
        assert (
            os.path.basename(ts1a_folder.crestopt_log_filepath)
            == "crestopt.log"
        )

        assert ts1a_folder.coord_filepath is not None
        assert os.path.basename(ts1a_folder.coord_filepath) == "coord"

    def test_folder_hexane_failed(self, crest_hexane_failed_outfolder):
        """Test CRESTFolder with a failed hexane conformational search."""
        assert os.path.exists(crest_hexane_failed_outfolder)
        hexane_folder = CRESTFolder(crest_hexane_failed_outfolder)
        assert hexane_folder.is_crest_calculation_directory

        assert hexane_folder.crest_out_filepath is not None
        assert (
            os.path.basename(hexane_folder.crest_out_filepath) == "hexane.out"
        )

        assert hexane_folder.conformers_xyz_filepath is None
        assert hexane_folder.best_xyz_filepath is None
        assert hexane_folder.rotamers_xyz_filepath is None
        assert hexane_folder.energies_filepath is None
        assert hexane_folder.constraints_inp_filepath is None

        assert hexane_folder.crestopt_log_filepath is not None
        assert (
            os.path.basename(hexane_folder.crestopt_log_filepath)
            == "crestopt.log"
        )

        assert hexane_folder.coord_filepath is not None
        assert os.path.basename(hexane_folder.coord_filepath) == "coord"


class TestCRESTOutput:
    """Tests for CRESTOutput class."""

    def test_octane_output(self, crest_octane_outfolder):
        """Test CRESTOutput class with octane conformational sampling output."""
        assert os.path.exists(crest_octane_outfolder)
        octane_output = CRESTOutput(folder=crest_octane_outfolder)
        assert octane_output.normal_termination
        assert len(octane_output.conformers) == 104
        assert octane_output.num_conformers == 104
        assert len(octane_output.rotamers) == 1342
        assert octane_output.num_rotamers == 1342

        assert octane_output.conformers[0].charge == 0
        assert octane_output.conformers[0].multiplicity == 1
        assert octane_output.best_conformer.charge == 0
        assert octane_output.best_conformer.multiplicity == 1
        assert octane_output.rotamers[0].charge == 0
        assert octane_output.rotamers[0].multiplicity == 1
        assert octane_output.conformers[0].num_atoms == 26
        assert octane_output.charge == 0
        assert octane_output.multiplicity == 1
        assert octane_output.num_atoms == 26

        assert np.isclose(octane_output.conformers[0].energy, -26.32263)
        assert np.isclose(octane_output.best_conformer.energy, -26.32263)
        assert np.isclose(octane_output.energies[0], -26.32263)
        assert np.isclose(octane_output.lowest_energy, -26.32263)

    def test_hexane_failed_output(self, crest_hexane_failed_outfolder):
        """Failed CREST run: no ensembles, no energies, abnormal termination."""
        assert os.path.exists(crest_hexane_failed_outfolder)
        hexane_output = CRESTOutput(folder=crest_hexane_failed_outfolder)

        assert not hexane_output.normal_termination
        assert hexane_output.conformers == []
        assert hexane_output.rotamers == []
        assert hexane_output.best_conformer is None
        assert hexane_output.energies == []
