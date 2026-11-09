"""Tests for reading real PySCF output files.

Every artifact under ``tests/data/PySCFTests/outputs`` was produced by
running ``chemsmart run pyscf`` against PySCF 2.14.0, so these tests check
the reader against bytes PySCF actually wrote rather than against output
from the writer that the reader is paired with.
"""

import os
import shutil

import h5py
import numpy as np
import pytest
from ase import units

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.utils.io import get_program_type_from_file
from tests.helpers import read_pyscf_output_folder


class TestPySCFArtifactProvenance:
    """The committed artifacts must be real calculations, not previews."""

    def test_artifacts_record_the_producing_pyscf_version(
        self,
        pyscf_water_sp_outfolder,
        pyscf_water_opt_outfolder,
        pyscf_water_hess_outfolder,
        pyscf_he_sp_outfolder,
    ):
        folders = (
            pyscf_water_sp_outfolder,
            pyscf_water_opt_outfolder,
            pyscf_water_hess_outfolder,
            pyscf_he_sp_outfolder,
        )
        for folder in folders:
            output = read_pyscf_output_folder(folder)
            assert output.version == "2.14.0"
            assert output.normal_termination
            assert output.engine == "cpu"

    def test_results_file_uses_navigable_top_level_groups(
        self, pyscf_water_opt_outfolder
    ):
        label = os.path.basename(pyscf_water_opt_outfolder)
        resultsfile = os.path.join(pyscf_water_opt_outfolder, f"{label}.h5")
        with h5py.File(resultsfile, "r") as handle:
            assert handle.attrs["schema_version"] == "2.0"
            for group in ("spec", "provenance", "status", "results"):
                assert group in handle


class TestPySCFSinglePointOutput:
    """Single point calculations on water in C-PCM water and on helium."""

    def test_water_single_point_energy_and_state(
        self, pyscf_water_sp_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_sp_outfolder)
        assert output.jobtype == "sp"
        assert output.method == "b3lyp"
        assert output.basis == "def2-svp"
        assert output.num_basis_functions == 24
        assert output.num_shells == 12
        assert output.charge == 0
        assert output.multiplicity == 1
        assert output.spin == 0
        assert output.chemical_symbols == ["O", "H", "H"]
        assert output.energies[-1] == -76.35923673423423

    def test_water_single_point_records_the_applied_solvent(
        self, pyscf_water_sp_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_sp_outfolder)
        assert output.solvent_on
        assert output.solvent_model == "cpcm"
        assert output.solvent_id == "water"
        assert output.freq is False
        assert output.vibrational_frequencies is None

    def test_water_single_point_population_analysis(
        self, pyscf_water_sp_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_sp_outfolder)
        charges = np.asarray(output.mulliken_atomic_charges)
        assert np.isclose(charges[0], -0.31530136, rtol=1e-6)
        assert np.allclose(charges[1], charges[2])
        # A neutral molecule's Mulliken charges must sum to zero.
        assert np.isclose(charges.sum(), 0.0, atol=1e-8)
        # Water is polar along the C2 axis only.
        dipole = np.asarray(output.dipole_moment)
        assert np.isclose(dipole[2], 2.04264644, rtol=1e-6)
        assert np.allclose(dipole[:2], 0.0, atol=1e-10)

    def test_helium_single_point_is_a_single_atom(self, pyscf_he_sp_outfolder):
        output = read_pyscf_output_folder(pyscf_he_sp_outfolder)
        assert output.chemical_symbols == ["He"]
        assert output.num_basis_functions == 5
        assert output.energies[-1] == -2.9071053229634183
        assert output.point_group == "SO3"
        assert output.rotational_symmetry_number == 1
        assert output.solvent_on is False
        assert output.vibrational_frequencies is None
        assert np.allclose(output.dipole_moment, 0.0, atol=1e-10)


class TestPySCFHessianOutput:
    """Hessian stages, both standalone and following an optimization."""

    def test_optimized_water_is_a_minimum(self, pyscf_water_opt_outfolder):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        assert output.jobtype == "opt"
        assert output.freq
        frequencies = np.asarray(output.vibrational_frequencies)
        assert len(frequencies) == 3
        assert np.all(frequencies > 0)
        assert np.allclose(
            frequencies,
            [1638.790195902117, 3791.255644245163, 3886.4637672678896],
            rtol=1e-6,
        )

    def test_optimization_lowers_the_energy(self, pyscf_water_opt_outfolder):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        energies = np.asarray(output.energies)
        assert len(energies) == 2
        assert energies[-1] < energies[0]
        assert energies[-1] == -76.35832577320055

    def test_water_symmetry_is_detected(self, pyscf_water_opt_outfolder):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        assert output.point_group == "C2v"
        assert output.rotational_symmetry_number == 2

    def test_standalone_hessian_uses_the_supplied_geometry(
        self, pyscf_water_hess_outfolder, pyscf_water_opt_outfolder
    ):
        hess = read_pyscf_output_folder(pyscf_water_hess_outfolder)
        opt = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        assert hess.jobtype == "hess"
        # The hess job never optimizes, so its energy is the optimization's
        # starting point, not its result.
        assert hess.energies[-1] == opt.energies[0]
        assert hess.energies[-1] > opt.energies[-1]

    def test_hessian_away_from_a_minimum_shifts_the_frequencies(
        self, pyscf_water_hess_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_hess_outfolder)
        frequencies = np.asarray(output.vibrational_frequencies)
        assert len(frequencies) == 3
        assert np.allclose(
            frequencies,
            [1217.7127703819842, 4559.4036539387835, 4732.154823050101],
            rtol=1e-6,
        )

    def test_hessian_is_symmetric_and_modes_are_stored(
        self, pyscf_water_opt_outfolder
    ):
        label = os.path.basename(pyscf_water_opt_outfolder)
        resultsfile = os.path.join(pyscf_water_opt_outfolder, f"{label}.h5")
        with h5py.File(resultsfile, "r") as handle:
            results = handle["results"]
            hessian = np.asarray(results["hessian"])
            reduced_masses = np.asarray(results["reduced_masses"])
            force_constants = np.asarray(results["force_constants"])
        assert hessian.shape == (3, 3, 3, 3)
        pairs = hessian.transpose(1, 0, 3, 2)
        assert np.allclose(hessian, pairs, rtol=2e-5, atol=1e-7)
        assert len(reduced_masses) == 3
        assert np.all(reduced_masses > 0)
        assert len(force_constants) == 3


class TestPySCFMoleculeRoundTrip:
    """Structures recovered from a results file."""

    def test_molecule_preserves_float64_geometry(
        self, pyscf_water_opt_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        molecule = output.get_molecule()
        assert molecule.num_atoms == 3
        assert molecule.empirical_formula == "H2O"
        assert molecule.charge == 0
        assert molecule.multiplicity == 1
        assert molecule.positions.dtype == np.float64
        assert np.allclose(molecule.positions, output.positions)

    def test_helium_molecule_round_trip(self, pyscf_he_sp_outfolder):
        output = read_pyscf_output_folder(pyscf_he_sp_outfolder)
        molecule = output.get_molecule()
        assert molecule.num_atoms == 1
        assert molecule.empirical_formula == "He"
        assert np.allclose(molecule.positions, 0.0)


class TestPySCFSymmetryFromDistortedStart:
    """Optimising from a distorted start must recover the true symmetry."""

    def test_near_symmetric_geometry_gets_the_true_symmetry_number(
        self, pyscf_water_opt_distorted_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_opt_distorted_outfolder)
        assert output.point_group == "C2v"
        assert output.rotational_symmetry_number == 2

    def test_free_energy_does_not_depend_on_the_starting_geometry(
        self, pyscf_water_opt_distorted_outfolder, pyscf_water_opt_outfolder
    ):
        distorted = read_pyscf_output_folder(
            pyscf_water_opt_distorted_outfolder
        )
        symmetric = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        # A wrong sigma would differ by RT ln 2 = 6.5e-4 Hartree.
        assert distorted.gibbs_free_energy == pytest.approx(
            symmetric.gibbs_free_energy, abs=1e-5
        )


class TestPySCFStructureAccess:
    """The stored structure behaves like other output readers' structures."""

    def test_optimized_structure_is_available_after_convergence(
        self, pyscf_water_opt_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        assert output.optimized_structure.is_optimized_structure
        assert len(output.all_structures) == 1

    def test_index_selects_among_stored_structures(
        self, pyscf_water_opt_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        assert isinstance(output.get_molecule(index=":"), list)
        assert output.get_molecule(index="1").num_atoms == 3
        with pytest.raises(IndexError):
            output.get_molecule(index="99")

    def test_dipole_and_point_group_reach_the_molecule(
        self, pyscf_water_opt_outfolder
    ):
        molecule = read_pyscf_output_folder(
            pyscf_water_opt_outfolder
        ).get_molecule()
        assert molecule.point_group == "C2V"
        assert molecule.dipole_moment_magnitude == pytest.approx(
            np.linalg.norm(molecule.dipole_moment)
        )


class TestPySCFUnrestrictedOutput:
    """UKS B3LYP/def2-SVP on the OH radical, computed with PySCF 2.14."""

    def test_spin_channels_are_separated(self, pyscf_oh_radical_outfolder):
        output = read_pyscf_output_folder(pyscf_oh_radical_outfolder)
        assert output.multiplicity == 2
        assert output.spin == 1
        assert output.results["mo_energy"].ndim == 2
        assert len(output.alpha_occ_eigenvalues) == 5
        assert len(output.beta_occ_eigenvalues) == 4
        assert output.num_unpaired_electrons == 1

    def test_orbital_energies_are_in_ev(self, pyscf_oh_radical_outfolder):
        output = read_pyscf_output_folder(pyscf_oh_radical_outfolder)
        assert output.alpha_homo_energy == pytest.approx(
            float(output.results["mo_energy"][0][4]) * units.Hartree
        )
        assert output.alpha_homo_energy != output.beta_homo_energy

    def test_linear_molecule_symmetry(self, pyscf_oh_radical_outfolder):
        output = read_pyscf_output_folder(pyscf_oh_radical_outfolder)
        assert output.point_group == "Coov"
        assert output.rotational_symmetry_number == 1


class TestPySCFUnconvergedScf:
    """A real SCF stopped after one cycle."""

    def test_does_not_read_as_complete(self, pyscf_water_scf_failed_outfolder):
        output = read_pyscf_output_folder(pyscf_water_scf_failed_outfolder)
        assert output.normal_termination is False
        assert output.status["stages"]["scf"]["converged"] is False

    def test_no_forces_from_an_unconverged_density(
        self, pyscf_water_scf_failed_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_scf_failed_outfolder)
        assert output.forces is None
        assert (
            output.status["skipped_properties"]["forces"]
            == "final SCF did not converge"
        )


class TestPySCFUnconvergedOptimisation:
    """A real optimisation stopped by its step limit, with freq requested."""

    def test_structure_is_not_labelled_optimized(
        self, pyscf_water_opt_unconverged_outfolder
    ):
        output = read_pyscf_output_folder(
            pyscf_water_opt_unconverged_outfolder
        )
        assert output.normal_termination is False
        assert output.optimizer_converged is False
        assert output.get_molecule().is_optimized_structure is False
        assert output.optimized_structure is None

    def test_hessian_was_skipped_so_no_thermochemistry_is_reported(
        self, pyscf_water_opt_unconverged_outfolder
    ):
        output = read_pyscf_output_folder(
            pyscf_water_opt_unconverged_outfolder
        )
        assert output.status["stages"]["hess"]["completed"] is False
        assert output.freq is False
        assert output.vibrational_frequencies is None
        assert output.gibbs_free_energy is None
        assert output.internal_energy is None
        assert output.thermal_enthalpy_correction is None


class TestPySCFRunWithoutGeometry:
    """A real SMD request on a PySCF build without that module."""

    def test_failure_is_recorded_with_its_stage(
        self, pyscf_water_smd_failed_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_smd_failed_outfolder)
        assert output.normal_termination is False
        assert output.failure["type"] == "RuntimeError"
        assert output.failure["stage"] == "scf"
        assert "SMD" in output.failure["message"]

    def test_requesting_a_structure_explains_the_failure(
        self, pyscf_water_smd_failed_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_smd_failed_outfolder)
        with pytest.raises(ValueError, match="SMD module is not available"):
            output.get_molecule()


class TestPySCFOutputDescriptions:
    def test_repr_names_the_calculation(self, pyscf_water_opt_outfolder):
        text = repr(read_pyscf_output_folder(pyscf_water_opt_outfolder))
        assert "jobtype='opt'" in text
        assert "normal_termination=True" in text

    def test_listed_indices_return_a_list(self, pyscf_water_opt_outfolder):
        selected = read_pyscf_output_folder(
            pyscf_water_opt_outfolder
        ).get_molecule(index="1,1")
        assert [molecule.num_atoms for molecule in selected] == [3, 3]

    def test_single_index_can_be_wrapped(self, pyscf_water_opt_outfolder):
        selected = read_pyscf_output_folder(
            pyscf_water_opt_outfolder
        ).get_molecule(index="-1", return_list=True)
        assert len(selected) == 1

    def test_single_point_reports_no_thermochemistry(
        self, pyscf_water_sp_outfolder
    ):
        output = read_pyscf_output_folder(pyscf_water_sp_outfolder)
        assert output.internal_energy is None
        assert output.thermal_gibbs_free_energy_correction is None

    def test_route_names_every_identity_setting(
        self, pyscf_water_sp_outfolder, tmp_path
    ):
        """No committed output used dispersion, as it needs an extra package."""
        label = os.path.basename(pyscf_water_sp_outfolder)
        for suffix in (".out", ".h5"):
            shutil.copy(
                os.path.join(pyscf_water_sp_outfolder, label + suffix),
                tmp_path / f"calc{suffix}",
            )
        with h5py.File(tmp_path / "calc.h5", "a") as handle:
            del handle["spec/dispersion"]
            handle["spec"].create_dataset("dispersion", data="d3bj")
        route = PySCFOutput(tmp_path / "calc.out").route_string
        assert "d3bj" in route
        assert "cpcm(water)" in route
        assert "engine=cpu" in route


class TestPySCFFailureReason:
    def test_recorded_exception_is_described(
        self, pyscf_water_smd_failed_outfolder
    ):
        reason = read_pyscf_output_folder(
            pyscf_water_smd_failed_outfolder
        ).failure_reason
        assert reason.startswith("RuntimeError at stage 'scf'")
        assert "SMD module is not available" in reason

    def test_unconverged_stages_are_named(
        self, pyscf_water_opt_unconverged_outfolder
    ):
        reason = read_pyscf_output_folder(
            pyscf_water_opt_unconverged_outfolder
        ).failure_reason
        assert "did not finish" in reason
        assert "opt" in reason
        assert "hess" in reason

    def test_successful_run_has_no_reason(self, pyscf_water_opt_outfolder):
        output = read_pyscf_output_folder(pyscf_water_opt_outfolder)
        assert output.failure_reason is None


class TestReadingResultsFilesDirectly:
    def _results_file(self, folder):
        return os.path.join(folder, f"{os.path.basename(folder)}.h5")

    def test_matches_reading_through_the_log(self, pyscf_water_opt_outfolder):
        log = self._results_file(pyscf_water_opt_outfolder)[:-3] + ".out"
        direct = Molecule.from_filepath(
            self._results_file(pyscf_water_opt_outfolder)
        )
        via_log = Molecule.from_filepath(log)
        np.testing.assert_array_equal(direct.positions, via_log.positions)
        assert direct.energy == via_log.energy

    def test_works_without_the_log(self, pyscf_water_opt_outfolder, tmp_path):
        path = tmp_path / "calc.h5"
        shutil.copy(self._results_file(pyscf_water_opt_outfolder), path)
        assert Molecule.from_filepath(str(path)).num_atoms == 3

    def test_index_and_return_list_are_honoured(
        self, pyscf_water_opt_outfolder
    ):
        path = self._results_file(pyscf_water_opt_outfolder)
        assert len(Molecule.from_filepath(path, index=":")) == 1
        assert len(Molecule.from_filepath(path, return_list=True)) == 1

    def test_failed_run_reports_its_failure(
        self, pyscf_water_smd_failed_outfolder
    ):
        with pytest.raises(ValueError, match="SMD module is not available"):
            Molecule.from_filepath(
                self._results_file(pyscf_water_smd_failed_outfolder)
            )

    def test_empty_log_reports_the_failure_beside_it(
        self, pyscf_water_smd_failed_outfolder, tmp_path
    ):
        shutil.copy(
            self._results_file(pyscf_water_smd_failed_outfolder),
            tmp_path / "calc.h5",
        )
        (tmp_path / "calc.out").write_text("")
        with pytest.raises(ValueError, match="SMD module is not available"):
            Molecule.from_filepath(str(tmp_path / "calc.out"))

    def test_empty_log_without_results_still_returns_none(self, tmp_path):
        (tmp_path / "calc.out").write_text("")
        assert Molecule.from_filepath(str(tmp_path / "calc.out")) is None

    def test_other_hdf5_files_are_refused_clearly(self, tmp_path):
        path = tmp_path / "other.h5"
        with h5py.File(path, "w") as handle:
            handle.create_dataset("data", data=[1.0, 2.0])
        with pytest.raises(ValueError, match="Unsupported PySCF results"):
            Molecule.from_filepath(str(path))

    def test_non_hdf5_file_is_refused_clearly(self, tmp_path):
        path = tmp_path / "other.h5"
        path.write_bytes(b"not an hdf5 file")
        with pytest.raises(ValueError, match="not a readable PySCF results"):
            Molecule.from_filepath(str(path))


class TestPySCFBenzeneOutput:
    """Optimisation and frequencies of benzene from a real PySCF run."""

    def test_run_finished_normally(self, pyscf_benzene_output):
        output = pyscf_benzene_output
        assert output.normal_termination
        assert output.optimizer_converged
        assert output.failure_reason is None
        assert output.spec["stages"] == ["scf", "opt", "hess"]
        assert (output.method, output.basis) == ("b3lyp", "def2-svp")
        assert output.version == "2.14.0"

    def test_optimisation_lowers_the_energy(self, pyscf_benzene_output):
        first, final = pyscf_benzene_output.energies
        assert final < first
        assert final == pytest.approx(-232.08478438, abs=1e-7)

    def test_geometry_is_a_regular_flat_hexagon(self, pyscf_benzene_output):
        positions = pyscf_benzene_output.positions
        symbols = pyscf_benzene_output.chemical_symbols
        distances = np.linalg.norm(
            positions[:, None] - positions[None], axis=-1
        )
        carbons = [i for i, s in enumerate(symbols) if s == "C"]
        hydrogens = [i for i, s in enumerate(symbols) if s == "H"]
        carbon_carbon = [
            distances[i, j]
            for i in carbons
            for j in carbons
            if i < j and distances[i, j] < 1.6
        ]
        carbon_hydrogen = [
            distances[i, j]
            for i in carbons
            for j in hydrogens
            if distances[i, j] < 1.3
        ]
        assert carbon_carbon == pytest.approx([1.3988] * 6, abs=1e-3)
        assert carbon_hydrogen == pytest.approx([1.0933] * 6, abs=1e-3)
        centred = positions - positions.mean(axis=0)
        assert np.linalg.svd(centred, compute_uv=False)[2] < 1e-3

    @pytest.mark.xfail(
        strict=True,
        reason="fixture was computed before the 1e-2 point-group tolerance",
    )
    def test_near_exact_symmetry_gives_the_true_symmetry_number(
        self, pyscf_benzene_output
    ):
        """C-C bonds within 2e-5 Angstrom once read as C2h, with sigma = 2."""
        assert pyscf_benzene_output.point_group == "D6h"
        assert pyscf_benzene_output.rotational_symmetry_number == 12

    @pytest.mark.xfail(
        strict=True,
        reason="fixture was computed before the 1e-2 point-group tolerance",
    )
    def test_free_energy_uses_the_true_symmetry_number(
        self, pyscf_benzene_output
    ):
        assert pyscf_benzene_output.gibbs_free_energy == pytest.approx(
            -232.009398, abs=1e-5
        )

    def test_minimum_has_all_real_frequencies(self, pyscf_benzene_output):
        frequencies = pyscf_benzene_output.vibrational_frequencies
        assert len(frequencies) == 30
        assert min(frequencies) == pytest.approx(418.5, abs=0.5)
        assert max(frequencies) == pytest.approx(3198.3, abs=0.5)

    def test_final_forces_vanish_at_the_minimum(self, pyscf_benzene_output):
        assert np.abs(pyscf_benzene_output.forces).max() < 1e-4
        assert pyscf_benzene_output.forces.shape == (12, 3)

    def test_nonpolar_molecule_has_no_dipole(self, pyscf_benzene_output):
        molecule = pyscf_benzene_output.get_molecule()
        assert molecule.dipole_moment_magnitude < 1e-3

    def test_frontier_orbitals_are_in_ev(self, pyscf_benzene_output):
        output = pyscf_benzene_output
        assert output.homo_energy == pytest.approx(-6.96678, abs=1e-4)
        assert output.lumo_energy == pytest.approx(-0.24861, abs=1e-4)
        assert output.fmo_gap == pytest.approx(6.71817, abs=1e-4)

    def test_zero_point_energy(self, pyscf_benzene_output):
        # Unlike the free energy, this does not depend on the point group.
        assert pyscf_benzene_output.zero_point_energy == pytest.approx(
            0.100499, abs=1e-5
        )

    def test_structure_reads_as_aromatic_benzene(self, pyscf_benzene_output):
        molecule = pyscf_benzene_output.get_molecule()
        assert molecule.chemical_formula == "C6H6"
        assert molecule.is_aromatic
        assert molecule.is_optimized_structure


class TestPySCFResultsWithoutTheirLog:
    def _only_results(self, folder, tmp_path):
        path = tmp_path / "calc.h5"
        shutil.copy(os.path.join(folder, "benzene_gas_phase.h5"), path)
        return path

    def test_thermochemistry_reads_the_results_file(
        self, pyscf_benzene_outfolder, pyscf_benzene_output, tmp_path
    ):
        path = self._only_results(pyscf_benzene_outfolder, tmp_path)
        assert PySCFOutput(path).gibbs_free_energy == pytest.approx(
            pyscf_benzene_output.gibbs_free_energy
        )

    def test_absent_log_does_not_stop_thermochemistry(
        self, pyscf_benzene_outfolder, pyscf_benzene_output, tmp_path
    ):
        self._only_results(pyscf_benzene_outfolder, tmp_path)
        output = PySCFOutput(tmp_path / "calc.out")
        assert output.gibbs_free_energy == pytest.approx(
            pyscf_benzene_output.gibbs_free_energy
        )

    def test_results_file_is_recognised_as_pyscf(
        self, pyscf_benzene_outfolder
    ):
        path = os.path.join(pyscf_benzene_outfolder, "benzene_gas_phase.h5")
        assert get_program_type_from_file(path) == "pyscf"

    def test_molecule_reads_the_results_file(self, pyscf_benzene_outfolder):
        path = os.path.join(pyscf_benzene_outfolder, "benzene_gas_phase.h5")
        assert Molecule.from_filepath(path).chemical_formula == "C6H6"


class TestPySCFLogRecognition:
    """A PySCF log is recognised whichever driver generation wrote it.

    The driver echo is longer than the 200-line detection window, so the
    driver's own signature line is the only thing that identifies the log.
    """

    def test_log_from_an_earlier_driver_is_recognised(
        self, pyscf_he_sp_outfolder
    ):
        path = os.path.join(pyscf_he_sp_outfolder, "he_sp_gas_phase.out")
        with open(path) as f:
            assert "# ChemSmart PySCF driver" in f.read(200)
        assert get_program_type_from_file(path) == "pyscf"
        molecule = Molecule.from_filepath(path)
        assert molecule.chemical_formula == "He"
        assert molecule.multiplicity == 1

    def test_log_from_the_current_driver_is_recognised(
        self, pyscf_benzene_outfolder
    ):
        path = os.path.join(pyscf_benzene_outfolder, "benzene_gas_phase.out")
        with open(path) as f:
            assert "# CHEMSMART PySCF driver" in f.read(200)
        assert get_program_type_from_file(path) == "pyscf"
