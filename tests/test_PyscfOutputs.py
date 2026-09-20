"""Tests for reading real PySCF output files.

Every artifact under ``tests/data/PySCFTests/outputs`` was produced by
running ``chemsmart run pyscf`` against PySCF 2.14.0, so these tests check
the reader against bytes PySCF actually wrote rather than against output
from the writer that the reader is paired with.
"""

import os

import h5py
import numpy as np

from chemsmart.io.pyscf.output import PySCFOutput


def _output(folder):
    """Return a PySCFOutput for the run stored in ``folder``."""
    label = os.path.basename(folder)
    logfile = os.path.join(folder, f"{label}.out")
    assert os.path.exists(logfile)
    assert os.path.exists(os.path.join(folder, f"{label}.h5"))
    return PySCFOutput(logfile)


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
            output = _output(folder)
            assert output.version == "2.14.0"
            assert output.normal_termination
            assert output.engine == "cpu"

    def test_results_file_uses_navigable_top_level_groups(
        self, pyscf_water_opt_outfolder
    ):
        label = os.path.basename(pyscf_water_opt_outfolder)
        resultsfile = os.path.join(pyscf_water_opt_outfolder, f"{label}.h5")
        with h5py.File(resultsfile, "r") as handle:
            assert handle.attrs["schema_version"] in ("1.0", "2.0")
            for group in ("spec", "provenance", "status", "results"):
                assert group in handle


class TestPySCFSinglePointOutput:
    """Single point calculations on water in C-PCM water and on helium."""

    def test_water_single_point_energy_and_state(
        self, pyscf_water_sp_outfolder
    ):
        output = _output(pyscf_water_sp_outfolder)
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
        output = _output(pyscf_water_sp_outfolder)
        assert output.solvent_on
        assert output.solvent_model == "cpcm"
        assert output.solvent_id == "water"
        assert output.freq is False
        assert output.vibrational_frequencies is None

    def test_water_single_point_population_analysis(
        self, pyscf_water_sp_outfolder
    ):
        output = _output(pyscf_water_sp_outfolder)
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
        output = _output(pyscf_he_sp_outfolder)
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
        output = _output(pyscf_water_opt_outfolder)
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
        output = _output(pyscf_water_opt_outfolder)
        energies = np.asarray(output.energies)
        assert len(energies) == 2
        assert energies[-1] < energies[0]
        assert energies[-1] == -76.35832577320055

    def test_water_symmetry_is_detected(self, pyscf_water_opt_outfolder):
        output = _output(pyscf_water_opt_outfolder)
        assert output.point_group == "C2v"
        assert output.rotational_symmetry_number == 2

    def test_standalone_hessian_uses_the_supplied_geometry(
        self, pyscf_water_hess_outfolder, pyscf_water_opt_outfolder
    ):
        hess = _output(pyscf_water_hess_outfolder)
        opt = _output(pyscf_water_opt_outfolder)
        assert hess.jobtype == "hess"
        # The hess job never optimizes, so its energy is the optimization's
        # starting point, not its result.
        assert hess.energies[-1] == opt.energies[0]
        assert hess.energies[-1] > opt.energies[-1]

    def test_hessian_away_from_a_minimum_shifts_the_frequencies(
        self, pyscf_water_hess_outfolder
    ):
        output = _output(pyscf_water_hess_outfolder)
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
        output = _output(pyscf_water_opt_outfolder)
        molecule = output.get_molecule()
        assert molecule.num_atoms == 3
        assert molecule.empirical_formula == "H2O"
        assert molecule.charge == 0
        assert molecule.multiplicity == 1
        assert molecule.positions.dtype == np.float64
        assert np.allclose(molecule.positions, output.positions)

    def test_helium_molecule_round_trip(self, pyscf_he_sp_outfolder):
        output = _output(pyscf_he_sp_outfolder)
        molecule = output.get_molecule()
        assert molecule.num_atoms == 1
        assert molecule.empirical_formula == "He"
        assert np.allclose(molecule.positions, 0.0)
