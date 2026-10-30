import os
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

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
        # CREGEN never ran for this failed job, so none of its summary
        # lines are present -- every one of these falls through to None.
        assert hexane_main_out.rmsd_threshold is None
        assert hexane_main_out.bconst_threshold is None
        assert hexane_main_out.population_threshold is None
        assert hexane_main_out.energy_window is None
        assert hexane_main_out.reference_state_energy is None
        assert hexane_main_out.temperature is None
        assert hexane_main_out.average_energy_ensemble is None
        assert hexane_main_out.entropy_ensemble is None
        assert hexane_main_out.free_energy_ensemble is None
        assert hexane_main_out.population_of_lowest is None

    def test_get_route_accepts_prompt_style_prefix(self, tmp_path):
        """Some CREST logs prefix the command line with '>' instead of
        '$' depending on how the run was invoked/piped."""
        out_file = tmp_path / "synthetic.out"
        out_file.write_text("some header\n> crest mol.xyz --gfn2\nmore text\n")
        main_out = CRESTMainOut(str(out_file))
        assert main_out.route_string == "crest mol.xyz --gfn2"

    def test_get_route_returns_none_when_no_command_line_present(
        self, tmp_path
    ):
        out_file = tmp_path / "synthetic.out"
        out_file.write_text("some header\nno command line here\n")
        main_out = CRESTMainOut(str(out_file))
        assert main_out.route_string is None

    def test_topology_mismatch_true_when_warning_present(self, tmp_path):
        out_file = tmp_path / "synthetic.out"
        out_file.write_text(
            "$ crest mol.xyz --gfn2\n"
            "*WARNING* Change in topology detected!\n"
        )
        main_out = CRESTMainOut(str(out_file))
        assert main_out.topology_mismatch is True

    def test_malformed_summary_values_fall_back_to_none(self, tmp_path):
        """Every CREGEN/ensemble-summary line present but with a
        non-numeric value hits the except-continue branch and, since no
        other matching line follows, ultimately returns None."""
        out_file = tmp_path / "synthetic.out"
        out_file.write_text(
            "$ crest mol.xyz --gfn2\n"
            " number of atoms                : notanumber\n"
            " RMSD threshold                 : notanumber\n"
            " Bconst threshold               : notanumber\n"
            " population threshold           : notanumber\n"
            " sorting energy window (EWIN)   : notanumber / kcal*mol-1\n"
            " reference state Etot           : notanumber\n"
            "T /K                            : notanumber\n"
            "E lowest                        :   notanumber\n"
            "ensemble average energy (kcal)  : notanumber\n"
            "ensemble entropy (J/mol K, cal/mol K) : notanumber notanumber\n"
            "ensemble free energy (kcal/mol) : notanumber\n"
            "population of lowest in %       : notanumber\n"
            " number of unique conformers for further calc   notanumber\n"
            " total number unique points considered further  notanumber\n"
        )
        main_out = CRESTMainOut(str(out_file))
        assert main_out.num_atoms is None
        assert main_out.rmsd_threshold is None
        assert main_out.bconst_threshold is None
        assert main_out.population_threshold is None
        assert main_out.energy_window is None
        assert main_out.reference_state_energy is None
        assert main_out.temperature is None
        assert main_out.lowest_energy is None
        # average/entropy/free energy all short-circuit to None as soon
        # as reference_state_energy is None, without even scanning for
        # their own markers -- covered by the reference_state_energy
        # malformed-value branch above already returning None.
        assert main_out.average_energy_ensemble is None
        assert main_out.entropy_ensemble is None
        assert main_out.free_energy_ensemble is None
        assert main_out.population_of_lowest is None
        assert main_out.num_conformers is None
        assert main_out.num_rotamers is None

    def test_entropy_and_average_energy_malformed_when_reference_present(
        self, tmp_path
    ):
        """With a valid reference_state_energy (so the early-return guard
        in average_energy_ensemble is skipped), a malformed ensemble
        average/entropy line still hits its own except-continue branch."""
        out_file = tmp_path / "synthetic.out"
        out_file.write_text(
            "$ crest mol.xyz --gfn2\n"
            " reference state Etot           : -1.0\n"
            "ensemble average energy (kcal)  : notanumber\n"
            "ensemble entropy (J/mol K, cal/mol K) : notanumber notanumber\n"
        )
        main_out = CRESTMainOut(str(out_file))
        assert main_out.reference_state_energy == -1.0
        assert main_out.average_energy_ensemble is None
        assert main_out.entropy_ensemble is None
        # free_energy_ensemble's own guard short-circuits to None too,
        # since it depends on average_energy_ensemble.
        assert main_out.free_energy_ensemble is None

    def test_free_energy_malformed_when_average_energy_present(self, tmp_path):
        """With a valid average_energy_ensemble (so free_energy_ensemble's
        own early-return guard is skipped), a malformed ensemble free
        energy line hits free_energy_ensemble's own except-continue
        branch instead of short-circuiting at the guard."""
        out_file = tmp_path / "synthetic.out"
        out_file.write_text(
            "$ crest mol.xyz --gfn2\n"
            " reference state Etot           : -1.0\n"
            "ensemble average energy (kcal)  : 1.0\n"
            "ensemble free energy (kcal/mol) : notanumber\n"
        )
        main_out = CRESTMainOut(str(out_file))
        assert main_out.average_energy_ensemble is not None
        assert main_out.free_energy_ensemble is None


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

    def test_skips_lines_with_fewer_than_two_columns(self, tmp_path):
        energies_file = tmp_path / "crest.energies"
        energies_file.write_text("1\n2 0.5\n")
        energies = CRESTEnergiesFile(str(energies_file))
        assert energies.relative_energies == [0.5]

    def test_skips_lines_with_non_numeric_energy(self, tmp_path):
        energies_file = tmp_path / "crest.energies"
        energies_file.write_text("1 notanumber\n2 0.5\n")
        energies = CRESTEnergiesFile(str(energies_file))
        assert energies.relative_energies == [0.5]

    def test_empty_file_gives_none_num_conformers_and_energy_range(
        self, tmp_path
    ):
        energies_file = tmp_path / "crest.energies"
        energies_file.write_text("")
        energies = CRESTEnergiesFile(str(energies_file))
        assert energies.relative_energies == []
        assert energies.num_conformers is None
        assert energies.energy_range is None


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

    def test_is_crest_calculation_directory_false_without_out_file(
        self, tmp_path
    ):
        folder = CRESTFolder(str(tmp_path))
        assert folder.crest_out_filepath is None
        assert folder.is_crest_calculation_directory is False

    def test_is_crest_calculation_directory_false_without_markers(
        self, tmp_path, crest_octane_outfolder
    ):
        """A folder with a valid CREST *.out file but none of the
        common auxiliary marker files is not considered a valid CREST
        calculation directory."""
        import shutil

        shutil.copy(
            os.path.join(crest_octane_outfolder, "octane_conformers.out"),
            tmp_path / "octane_conformers.out",
        )
        folder = CRESTFolder(str(tmp_path))
        assert folder.crest_out_filepath is not None
        assert folder.is_crest_calculation_directory is False

    def test_crest_out_filepath_raises_on_multiple_output_files(
        self, tmp_path, crest_octane_outfolder
    ):
        import shutil

        src = os.path.join(crest_octane_outfolder, "octane_conformers.out")
        shutil.copy(src, tmp_path / "octane_conformers.out")
        shutil.copy(src, tmp_path / "octane_conformers_copy.out")
        folder = CRESTFolder(str(tmp_path))
        with pytest.raises(ValueError, match="Multiple CREST main output"):
            folder.crest_out_filepath


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

    @staticmethod
    def _make_output_with_mock_folder(**folder_overrides):
        output = CRESTOutput.__new__(CRESTOutput)
        mock_folder = MagicMock()
        for attr in (
            "crest_out_filepath",
            "energies_filepath",
            "conformers_xyz_filepath",
            "rotamers_xyz_filepath",
            "best_xyz_filepath",
        ):
            setattr(mock_folder, attr, folder_overrides.get(attr))
        output.folder = mock_folder
        return output

    def test_normal_termination_false_without_main_out(self):
        output = self._make_output_with_mock_folder()
        assert output.main_out is None
        assert output.normal_termination is False

    def test_conformers_empty_list_on_read_exception(self):
        output = self._make_output_with_mock_folder(
            conformers_xyz_filepath="/some/conformers.xyz"
        )
        with patch(
            "chemsmart.io.crest.output.Molecule.from_filepath",
            side_effect=ValueError("bad xyz"),
        ):
            assert output.conformers == []

    def test_rotamers_empty_list_on_read_exception(self):
        output = self._make_output_with_mock_folder(
            rotamers_xyz_filepath="/some/rotamers.xyz"
        )
        with patch(
            "chemsmart.io.crest.output.Molecule.from_filepath",
            side_effect=ValueError("bad xyz"),
        ):
            assert output.rotamers == []

    def test_best_conformer_falls_back_to_first_conformer_when_missing(
        self, crest_octane_outfolder
    ):
        output = self._make_output_with_mock_folder(
            crest_out_filepath=os.path.join(
                crest_octane_outfolder, "octane_conformers.out"
            ),
            conformers_xyz_filepath="/some/conformers.xyz",
        )
        fake_molecule = MagicMock()
        with patch(
            "chemsmart.io.crest.output.Molecule.from_filepath",
            return_value=[fake_molecule],
        ):
            assert output.best_conformer is fake_molecule

    def test_best_conformer_none_when_no_best_file_and_no_conformers(self):
        output = self._make_output_with_mock_folder()
        assert output.best_conformer is None

    def test_best_conformer_none_on_read_exception(self):
        output = self._make_output_with_mock_folder(
            best_xyz_filepath="/some/best.xyz"
        )
        with patch(
            "chemsmart.io.crest.output.Molecule.from_filepath",
            side_effect=ValueError("bad xyz"),
        ):
            assert output.best_conformer is None

    def test_energies_falls_back_to_energies_file_and_main_out(
        self, crest_octane_outfolder
    ):
        """When conformers carry no per-molecule energy, fall back to
        main_out.lowest_energy + energies_file.relative_energies."""
        output = self._make_output_with_mock_folder(
            crest_out_filepath=os.path.join(
                crest_octane_outfolder, "octane_conformers.out"
            ),
            energies_filepath=os.path.join(
                crest_octane_outfolder, "crest.energies"
            ),
        )
        # conformers_xyz_filepath is None -> self.conformers == []
        assert output.conformers == []
        assert output.main_out is not None
        assert output.energies_file is not None
        assert output.energies
        assert np.isclose(output.energies[0], output.main_out.lowest_energy)

    def test_num_conformers_falls_back_to_energies_file(
        self, crest_octane_outfolder
    ):
        output = self._make_output_with_mock_folder(
            energies_filepath=os.path.join(
                crest_octane_outfolder, "crest.energies"
            ),
        )
        assert output.main_out is None
        assert output.num_conformers == output.energies_file.num_conformers

    def test_num_conformers_falls_back_to_len_conformers(
        self, crest_hexane_failed_outfolder
    ):
        # hexane's main_out has no "number of unique conformers" line, so
        # main_out.num_conformers is None and no energies_file exists
        # either -- forces the len(self.conformers) fallback.
        output = self._make_output_with_mock_folder(
            crest_out_filepath=os.path.join(
                crest_hexane_failed_outfolder, "hexane.out"
            ),
            conformers_xyz_filepath="/some/conformers.xyz",
        )
        assert output.main_out.num_conformers is None
        fake_molecule = MagicMock()
        with patch(
            "chemsmart.io.crest.output.Molecule.from_filepath",
            return_value=[fake_molecule],
        ):
            assert output.num_conformers == 1

    def test_num_rotamers_falls_back_to_len_rotamers(
        self, crest_hexane_failed_outfolder
    ):
        # hexane's main_out has no "total number unique points" line, so
        # main_out.num_rotamers is None -- forces the len(self.rotamers)
        # fallback (no rotamers_xyz_filepath here, so it's 0).
        output = self._make_output_with_mock_folder(
            crest_out_filepath=os.path.join(
                crest_hexane_failed_outfolder, "hexane.out"
            ),
        )
        assert output.main_out.num_rotamers is None
        assert output.num_rotamers == 0


class TestCRESTRouteInit:
    def test_lowercases_and_splits_route_string(self):
        route = CRESTRoute("CREST Mol.xyz --GFN2 --Chrg 0")
        assert route.route_string == "crest mol.xyz --gfn2 --chrg 0"
        assert route.route_inputs == [
            "crest",
            "mol.xyz",
            "--gfn2",
            "--chrg",
            "0",
        ]


class TestCRESTRouteGfnVersionAndMethod:
    @pytest.mark.parametrize(
        "flag,expected",
        [
            ("--gfn1", "gfn1"),
            ("--gfn2", "gfn2"),
            ("--gfnff", "gfnff"),
            ("--gfn2//gfnff", "gfn2//gfnff"),
        ],
    )
    def test_explicit_gfn_flags(self, flag, expected):
        route = CRESTRoute(f"crest mol.xyz {flag}")
        assert route.gfn_version == expected
        assert route.method == expected

    def test_gff_alias_for_gfnff(self):
        route = CRESTRoute("crest mol.xyz --gff")
        assert route.gfn_version == "gfnff"

    @pytest.mark.parametrize(
        "num,expected",
        [("0", "gfn0"), ("1", "gfn1"), ("2", "gfn2")],
    )
    def test_gfn_int_format(self, num, expected):
        route = CRESTRoute(f"crest mol.xyz --gfn {num}")
        assert route.gfn_version == expected

    def test_gfn_int_out_of_supported_range_returns_none(self):
        route = CRESTRoute("crest mol.xyz --gfn 3")
        assert route.gfn_version is None

    def test_gfn_non_integer_value_returns_none(self):
        route = CRESTRoute("crest mol.xyz --gfn abc")
        assert route.gfn_version is None

    def test_gfn_flag_at_end_of_route_returns_none(self):
        route = CRESTRoute("crest mol.xyz --gfn")
        assert route.gfn_version is None

    def test_no_method_specified_returns_none(self):
        route = CRESTRoute("crest mol.xyz")
        assert route.gfn_version is None
        assert route.method is None


class TestCRESTRouteBasis:
    def test_basis_always_default(self):
        assert CRESTRoute("crest mol.xyz").basis == "default"


class TestCRESTRouteOptimizationLevel:
    @pytest.mark.parametrize(
        "flag",
        ["--optlev", "-optlev", "--opt", "-o", "--optlevel", "--ancopt"],
    )
    def test_each_flag_reads_following_level(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} tight")
        assert route.optimization_level == "tight"

    def test_verytight_normalized_to_vtight(self):
        route = CRESTRoute("crest mol.xyz --opt verytight")
        assert route.optimization_level == "vtight"

    def test_unknown_level_returns_none(self):
        route = CRESTRoute("crest mol.xyz --opt bogus_level")
        assert route.optimization_level is None

    def test_flag_at_end_of_route_returns_none(self):
        route = CRESTRoute("crest mol.xyz --opt")
        assert route.optimization_level is None

    def test_no_opt_flag_returns_none(self):
        assert CRESTRoute("crest mol.xyz").optimization_level is None


class TestCRESTRouteChargeAndUhf:
    @pytest.mark.parametrize("flag", ["--chrg", "-chrg"])
    def test_charge_flags(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} -1")
        assert route.charge == -1

    def test_charge_invalid_value_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --chrg notanumber")
        assert route.charge == 0

    def test_charge_flag_at_end_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --chrg")
        assert route.charge == 0

    def test_no_charge_flag_defaults_to_zero(self):
        assert CRESTRoute("crest mol.xyz").charge == 0

    @pytest.mark.parametrize("flag", ["--uhf", "-uhf"])
    def test_uhf_flags(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} 2")
        assert route.uhf == 2
        assert route.multiplicity == 3

    def test_uhf_invalid_value_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --uhf notanumber")
        assert route.uhf == 0

    def test_uhf_flag_at_end_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz --uhf")
        assert route.uhf == 0

    def test_no_uhf_flag_defaults_to_zero(self):
        route = CRESTRoute("crest mol.xyz")
        assert route.uhf == 0
        assert route.multiplicity == 1


class TestCRESTRouteJobtype:
    def test_jobtype_is_always_conformers(self):
        assert CRESTRoute("crest mol.xyz").jobtype == "conformers"
        assert CRESTRoute("crest mol.xyz --gfn2").get_jobtype() == "conformers"


class TestCRESTRouteSolventModelAndId:
    @pytest.mark.parametrize("model", ["gbsa", "alpb"])
    def test_each_solvent_model_flag(self, model):
        route = CRESTRoute(f"crest mol.xyz --{model} water")
        assert route.solvent_model == model
        assert route.get_solvent_model() == model

    def test_g_alias_for_gbsa(self):
        route = CRESTRoute("crest mol.xyz -g water")
        assert route.solvent_model == "gbsa"

    def test_no_solvent_model_returns_none(self):
        assert CRESTRoute("crest mol.xyz").solvent_model is None

    def test_solvent_id_extracted_when_model_present(self):
        route = CRESTRoute("crest mol.xyz --alpb toluene")
        assert route.solvent_id == "toluene"

    def test_solvent_id_none_when_no_solvent_model(self):
        route = CRESTRoute("crest mol.xyz toluene")
        assert route.solvent_id is None

    def test_solvent_id_none_when_model_present_but_no_matching_token(self):
        route = CRESTRoute("crest mol.xyz --alpb")
        assert route.solvent_id is None


class TestCRESTRouteBooleanFlags:
    @pytest.mark.parametrize(
        "attr,long_flag,short_flag",
        [
            ("nci", "--nci", "-nci"),
            ("constrained", "--cinp", "-cinp"),
            (
                "no_reference_topology_check",
                "--noreftopo",
                "-noreftopo",
            ),
            ("no_topology_check", "--notopo", "-notopo"),
        ],
    )
    def test_flag_true_for_long_and_short_forms(
        self, attr, long_flag, short_flag
    ):
        assert getattr(CRESTRoute(f"crest mol.xyz {long_flag}"), attr) is True
        assert getattr(CRESTRoute(f"crest mol.xyz {short_flag}"), attr) is True

    @pytest.mark.parametrize(
        "attr",
        [
            "nci",
            "constrained",
            "no_reference_topology_check",
            "no_topology_check",
        ],
    )
    def test_flag_false_when_absent(self, attr):
        assert getattr(CRESTRoute("crest mol.xyz"), attr) is False


class TestCRESTRouteThreads:
    @pytest.mark.parametrize("flag", ["-T", "--T"])
    def test_threads_flags_case_insensitive(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} 8")
        assert route.threads == 8

    def test_threads_invalid_value_returns_none(self):
        route = CRESTRoute("crest mol.xyz -T notanumber")
        assert route.threads is None

    def test_threads_flag_at_end_returns_none(self):
        route = CRESTRoute("crest mol.xyz -T")
        assert route.threads is None

    def test_no_threads_flag_returns_none(self):
        assert CRESTRoute("crest mol.xyz").threads is None


class TestCRESTRouteQuickMode:
    @pytest.mark.parametrize("mode", ["quick", "squick", "mquick"])
    @pytest.mark.parametrize("prefix", ["--", "-"])
    def test_quick_mode_flags(self, mode, prefix):
        route = CRESTRoute(f"crest mol.xyz {prefix}{mode}")
        assert route.quick_mode == mode

    def test_no_quick_mode_returns_none(self):
        assert CRESTRoute("crest mol.xyz").quick_mode is None


class TestCRESTRouteNumericThresholds:
    """Every 'flag REAL/INT' property sharing the same
    present/invalid/flag-at-end/absent shape."""

    _CASES = [
        ("rmsd_threshold", "--rthr", "-rthr", "0.5", 0.5),
        ("energy_threshold", "--ethr", "-ethr", "0.1", 0.1),
        ("bconst_threshold", "--bthr", "-bthr", "0.02", 0.02),
        ("population_threshold", "--pthr", "-pthr", "0.05", 0.05),
        ("temperature", "--temp", "-temp", "298.15", 298.15),
        ("shake", "--shake", "-shake", "1", 1),
        ("md_timestep", "--tstep", "-tstep", "4", 4),
        ("md_dump_step", "--mddump", "-mddump", "50", 50),
        ("vbias_dump_interval", "--vbdump", "-vbdump", "1.0", 1.0),
        ("additional_md_temperature", "--tnmd", "-tnmd", "400", 400.0),
        ("energy_window", "--ewin", "-ewin", "6.0", 6.0),
    ]

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_long_flag_parses_value(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {long_flag} {raw}")
        assert getattr(route, attr) == pytest.approx(expected)

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_short_flag_parses_value(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {short_flag} {raw}")
        assert getattr(route, attr) == pytest.approx(expected)

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_invalid_value_returns_none(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {long_flag} notanumber")
        assert getattr(route, attr) is None

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_flag_at_end_returns_none(
        self, attr, long_flag, short_flag, raw, expected
    ):
        route = CRESTRoute(f"crest mol.xyz {long_flag}")
        assert getattr(route, attr) is None

    @pytest.mark.parametrize("attr,long_flag,short_flag,raw,expected", _CASES)
    def test_absent_returns_none(
        self, attr, long_flag, short_flag, raw, expected
    ):
        assert getattr(CRESTRoute("crest mol.xyz"), attr) is None


class TestCRESTRouteMdLength:
    @pytest.mark.parametrize("flag", ["--mdlen", "-mdlen", "--len", "-len"])
    def test_numeric_value_parsed_as_float(self, flag):
        route = CRESTRoute(f"crest mol.xyz {flag} 50")
        assert route.md_length == pytest.approx(50.0)

    def test_x_prefixed_multiplicative_factor_returned_as_string(self):
        route = CRESTRoute("crest mol.xyz --mdlen x2.0")
        assert route.md_length == "x2.0"

    def test_non_numeric_non_x_value_returned_as_string(self):
        route = CRESTRoute("crest mol.xyz --mdlen notanumber")
        assert route.md_length == "notanumber"

    def test_flag_at_end_returns_none(self):
        assert CRESTRoute("crest mol.xyz --mdlen").md_length is None

    def test_absent_returns_none(self):
        assert CRESTRoute("crest mol.xyz").md_length is None
