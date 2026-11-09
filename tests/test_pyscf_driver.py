"""Run the generated PySCF driver for real, on molecules small enough to be fast.

Skipped where PySCF is not installed, such as Windows, which it does not support.
Every number here comes from an actual calculation.
"""

import os
from unittest.mock import patch

import numpy as np
import pytest

pytest.importorskip("pyscf")

from pyscf import __config__, cc, fci, gto, mp, symm  # noqa: E402
from pyscf.scf import hf  # noqa: E402

from chemsmart.database.assemble import PySCFAssembler  # noqa: E402
from chemsmart.io.gaussian.output import Gaussian16Output  # noqa: E402
from chemsmart.io.molecules.structure import Molecule  # noqa: E402
from chemsmart.io.orca.output import ORCAOutput  # noqa: E402
from chemsmart.io.pyscf import driver  # noqa: E402
from chemsmart.jobs.pyscf.hess import PySCFHessJob  # noqa: E402
from chemsmart.jobs.pyscf.irc import PySCFIRCJob  # noqa: E402
from chemsmart.jobs.pyscf.opt import PySCFOptJob  # noqa: E402
from chemsmart.jobs.pyscf.settings import PySCFJobSettings  # noqa: E402
from chemsmart.jobs.pyscf.singlepoint import (  # noqa: E402
    PySCFSinglePointJob,
)
from chemsmart.jobs.pyscf.ts import PySCFTSJob  # noqa: E402
from tests.helpers import run_pyscf_ccpvdz_sp  # noqa: E402

HF_MINIMAL = {"ab_initio": "hf", "basis": "sto-3g", "density_fit": False}
DFT_MINIMAL = {"basis": "sto-3g", "density_fit": False}
CORRELATED = {"density_fit": False, "scf_tol": 1e-11}


class TestSinglePoint:
    def test_energy_and_reader_agree_with_pyscf(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob, pyscf_h2_molecule, "sp", **HF_MINIMAL
        )
        assert returncode == 0
        assert output.normal_termination
        assert output.energies == pytest.approx([-1.116759], abs=1e-5)
        assert output.point_group == "Dooh"
        assert output.rotational_symmetry_number == 2
        assert output.spec["num_electrons"] == 2

    def test_single_point_computes_no_gradient_by_default(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob, pyscf_h2_molecule, "sp", **HF_MINIMAL
        )
        assert output.forces is None

    def test_requested_forces_have_the_physical_sign(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        # 0.74 Angstrom is longer than the HF/STO-3G bond (about 0.71), so a
        # force, not a gradient, pulls the second atom back towards the first.
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            forces=True,
            **HF_MINIMAL,
        )
        assert output.forces[1][2] < 0
        assert output.forces.sum(axis=0) == pytest.approx(
            np.zeros(3), abs=1e-6
        )

    def test_open_shell_molecule_uses_an_unrestricted_reference(
        self, run_pyscf_driver, pyscf_oh_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_oh_molecule,
            "sp",
            multiplicity=2,
            **HF_MINIMAL,
        )
        assert output.normal_termination
        assert output.results["mo_energy"].ndim == 2
        assert output.spec["spin"] == 1
        assert output.spec["nelec"] == [5, 4]
        assert len(output.alpha_occ_eigenvalues) == 5
        assert len(output.beta_occ_eigenvalues) == 4

    def test_dft_with_density_fitting_and_a_named_grid(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional="pbe",
            basis="sto-3g",
            density_fit=True,
            defgrid="defgrid1",
            scf_tol=1e-8,
            scf_maxiter=50,
        )
        assert output.normal_termination
        assert output.spec["atom_grid"] == [50, 194]
        assert output.spec["xc"] == "pbe"


class TestSolvent:
    def test_pcm_applies_the_requested_dielectric(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        _, gas = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            label="gas",
            **HF_MINIMAL,
        )
        _, toluene = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            label="toluene",
            solvent_model="cpcm",
            solvent_id="toluene",
            **HF_MINIMAL,
        )
        assert toluene.spec["solvent_eps"] == pytest.approx(2.3741)
        assert toluene.spec["solvent_method"] == "C-PCM"
        # Solvation stabilises a polar molecule.
        assert toluene.energies[-1] < gas.energies[-1]

    def test_unknown_solvent_is_recorded_as_a_failure(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            solvent_model="cpcm",
            solvent_id="not_a_solvent",
            **HF_MINIMAL,
        )
        assert returncode == 1
        assert not output.normal_termination
        assert "not_a_solvent" in output.failure["message"]

    def test_smd_runs_or_reports_the_missing_build_option(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        """PyPI wheels omit SMD; either outcome must be recorded faithfully."""
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            solvent_model="smd",
            solvent_id="water",
            **HF_MINIMAL,
        )
        assert output.spec["solvent_eps"] == pytest.approx(78.3553, abs=0.01)
        if returncode == 0:
            assert output.normal_termination
        else:
            assert "SMD" in output.failure["message"]


def test_dispersion_runs_or_reports_the_missing_package(
    run_pyscf_driver, pyscf_h2_molecule
):
    returncode, output = run_pyscf_driver(
        PySCFSinglePointJob,
        pyscf_h2_molecule,
        "sp",
        functional="pbe",
        basis="sto-3g",
        dispersion="d3bj",
    )
    if returncode == 0:
        assert output.normal_termination
    else:
        assert output.failure["type"] in ("RuntimeError", "TypeError")


class TestFailureHandling:
    def test_unconverged_scf_stores_no_forces_and_skips_the_hessian(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        returncode, output = run_pyscf_driver(
            PySCFHessJob,
            pyscf_water_molecule,
            "hess",
            scf_maxiter=1,
            forces=True,
            **HF_MINIMAL,
        )
        assert returncode == 0
        assert not output.normal_termination
        assert output.status["stages"]["scf"]["converged"] is False
        assert output.status["stages"]["hess"]["completed"] is False
        assert output.vibrational_frequencies is None
        assert output.forces is None
        assert "forces" in output.status["skipped_properties"]

    @pytest.mark.parametrize(
        ("override", "message"),
        [
            ({"stages": ["scf", "bogus"]}, "Unknown stage"),
            (
                {"stages": ["opt"], "opt_solver": "bogus", "opt_maxsteps": 1},
                "Unknown opt_solver",
            ),
        ],
    )
    def test_unknown_requests_are_recorded_not_raised(
        self, run_pyscf_driver, pyscf_h2_molecule, override, message
    ):
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            config_overrides=override,
            **HF_MINIMAL,
        )
        assert returncode == 1
        assert message in output.failure["message"]


class TestOptimisation:
    def test_geometric_lowers_the_energy_and_converges(
        self, run_pyscf_driver, pyscf_distorted_water_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFOptJob,
            pyscf_distorted_water_molecule,
            "opt",
            freq=False,
            **HF_MINIMAL,
        )
        assert output.normal_termination
        assert output.optimizer_converged
        assert output.energies[-1] < output.energies[0]
        assert output.get_molecule().is_optimized_structure
        # The stored forces belong to the converged geometry.
        assert np.abs(output.forces).max() < 1e-2

    def test_step_limit_is_not_reported_as_an_optimised_structure(
        self, run_pyscf_driver, pyscf_distorted_water_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFOptJob,
            pyscf_distorted_water_molecule,
            "opt",
            freq=True,
            opt_maxsteps=1,
            **HF_MINIMAL,
        )
        assert not output.normal_termination
        assert output.optimized_structure is None
        assert output.status["stages"]["hess"]["completed"] is False

    def test_optimisation_followed_by_frequencies(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFOptJob, pyscf_h2_molecule, "opt", freq=True, **HF_MINIMAL
        )
        assert output.normal_termination
        assert output.freq
        assert len(output.vibrational_frequencies) == 1
        assert output.vibrational_frequencies[0] > 0
        assert output.gibbs_free_energy is not None

    def test_ase_solver(self, run_pyscf_driver, pyscf_h2_molecule):
        _, output = run_pyscf_driver(
            PySCFOptJob,
            pyscf_h2_molecule,
            "opt",
            freq=False,
            opt_solver="ase",
            **HF_MINIMAL,
        )
        assert output.normal_termination
        assert output.spec["opt_solver"] == "ase"

    def test_berny_solver(self, run_pyscf_driver, pyscf_h2_molecule):
        pytest.importorskip("berny")
        _, output = run_pyscf_driver(
            PySCFOptJob,
            pyscf_h2_molecule,
            "opt",
            freq=False,
            opt_solver="berny",
            **HF_MINIMAL,
        )
        assert output.normal_termination


def test_imaginary_modes_are_stored_as_negative_wavenumbers(
    run_pyscf_driver, pyscf_linear_water_molecule
):
    """Linear water is a saddle point, so its bend must come out negative."""
    _, output = run_pyscf_driver(
        PySCFHessJob, pyscf_linear_water_molecule, "hess", **HF_MINIMAL
    )
    assert output.normal_termination
    assert min(output.vibrational_frequencies) < 0
    assert output.results["vibrational_frequencies"].dtype.kind == "f"


class TestPointGroupDetection:
    def _near_symmetric_water(self):
        return gto.M(
            atom=[
                ("O", (0.0, 0.0, 0.1173)),
                ("H", (0.0, 0.7572, -0.4692)),
                ("H", (0.0, -0.7573, -0.4692)),
            ],
            basis="sto-3g",
            verbose=0,
        )

    def test_tolerance_recovers_symmetry_of_a_converged_geometry(self):
        """A 1e-4 Angstrom asymmetry is typical of a converged optimisation."""
        mol = self._near_symmetric_water()
        assert symm.detect_symm(mol._atom)[0] == "Cs"
        assert driver._detect_point_group(mol, 1e-3) == "C2v"

    def test_pyscf_default_tolerance_is_restored(self):
        default = symm.geom.TOLERANCE
        driver._detect_point_group(self._near_symmetric_water(), 1e-3)
        assert symm.geom.TOLERANCE == default


class TestDefaultPointGroupTolerance:
    """The default must recover the symmetry of real optimised geometries."""

    def _default(self):
        return PySCFJobSettings().point_group_tol

    def test_benzene_from_a_real_optimisation(self, pyscf_benzene_output):
        mol = gto.M(
            atom=list(
                zip(
                    pyscf_benzene_output.chemical_symbols,
                    map(tuple, pyscf_benzene_output.positions),
                )
            ),
            basis="sto-3g",
            verbose=0,
        )
        assert driver._detect_point_group(mol, self._default()) == "D6h"
        # The previous default read this geometry as C2h, so sigma was 2.
        assert driver._detect_point_group(mol, 1e-3) == "C2h"

    @pytest.mark.parametrize(
        ("symbols", "positions", "expected"),
        [
            (
                ["C", "H", "H", "H", "H"],
                [
                    [0.0, 0.0, 0.0],
                    [0.66, 0.60, 0.62],
                    [-0.60, -0.66, 0.63],
                    [-0.62, 0.63, -0.60],
                    [0.63, -0.62, -0.66],
                ],
                "Td",
            ),
            (
                ["O", "C", "O"],
                [[0.0, 0.0, -1.17], [0.0, 0.01, 0.0], [0.0, 0.0, 1.18]],
                "Dooh",
            ),
        ],
    )
    def test_small_molecules_optimised_from_distorted_starts(
        self, run_pyscf_driver, symbols, positions, expected
    ):
        molecule = Molecule(
            symbols=symbols,
            positions=np.array(positions),
            charge=0,
            multiplicity=1,
        )
        _, output = run_pyscf_driver(
            PySCFOptJob, molecule, "opt", freq=False, **HF_MINIMAL
        )
        assert output.optimizer_converged
        assert output.point_group == expected

    def test_real_distortion_is_not_hidden(self):
        """Alternating ring radii of 0.02 Angstrom are a genuine D3h ring."""
        atoms = []
        for i in range(6):
            angle = i * np.pi / 3
            radius = 1.3988 + 0.01 * (1 if i % 2 else -1)
            atoms.append(
                ("C", (radius * np.cos(angle), radius * np.sin(angle), 0.0))
            )
            atoms.append(
                ("H", (2.4921 * np.cos(angle), 2.4921 * np.sin(angle), 0.0))
            )
        mol = gto.M(atom=atoms, basis="sto-3g", verbose=0)
        assert driver._detect_point_group(mol, self._default()) == "D3h"


def test_gpu_conversion_happens_before_the_solvent_is_attached(
    pyscf_h2_molecule, pyscf_driver_config
):
    """No GPU in CI: ``to_gpu`` is replaced, but the rest is real PySCF.

    gpu4pyscf needs this order, otherwise a CPU solvent stays on a GPU method.
    """
    attached_at_conversion = []

    def fake_to_gpu(method):
        attached_at_conversion.append(hasattr(method, "with_solvent"))
        return method

    config = pyscf_driver_config(
        PySCFSinglePointJob,
        pyscf_h2_molecule,
        "sp",
        engine="gpu",
        solvent_model="cpcm",
        solvent_id="water",
        **HF_MINIMAL,
    )
    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", verbose=0)
    with patch.object(hf.RHF, "to_gpu", fake_to_gpu):
        method = driver._build_method(config, mol)
    assert attached_at_conversion == [False]
    assert method.with_solvent.method == "C-PCM"


class TestAppliedMethodIsRecorded:
    """The results file states what PySCF computed, not only what was asked."""

    def test_b3lyp_is_the_gaussian_vwn_rpa_functional(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional="b3lyp",
            **DFT_MINIMAL,
        )
        applied = output.xc_applied
        assert output.reference == "RKS"
        assert applied["backend"] == "libxc"
        # LibXC 402 is HYB_GGA_XC_B3LYP, built on VWN_RPA as in Gaussian.
        assert applied["component_ids"] == [402]
        assert applied["hf_exchange_short_range"] == pytest.approx(0.2)
        assert applied["hf_exchange_long_range"] == pytest.approx(0.2)
        assert applied["range_separation_omega"] == 0.0
        assert applied["nonlocal_correlation"] is False

    def test_b3lyp5_is_composed_with_vwn5(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional="b3lyp5",
            **DFT_MINIMAL,
        )
        applied = output.xc_applied
        # Slater, B88, LYP and VWN5 (LibXC 1, 106, 131, 7), as in ORCA.
        assert applied["component_ids"] == [1, 106, 131, 7]
        assert applied["component_coefficients"] == pytest.approx(
            [0.08, 0.72, 0.81, 0.19]
        )
        assert applied["hf_exchange_short_range"] == pytest.approx(0.2)

    def test_range_separated_hybrid_records_both_exchange_ranges(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional="cam-b3lyp",
            **DFT_MINIMAL,
        )
        applied = output.xc_applied
        assert applied["range_separation_omega"] == pytest.approx(0.33)
        assert applied["hf_exchange_short_range"] == pytest.approx(0.19)
        assert applied["hf_exchange_long_range"] == pytest.approx(0.65)

    def test_defaults_left_unset_are_recorded_as_applied(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional="b3lyp",
            basis="def2-svp",
            defgrid="defgrid2",
        )
        numerics = output.numerics_applied
        assert output.spec["scf_tol"] is None
        assert numerics["conv_tol"] == pytest.approx(1e-9)
        assert numerics["grid"]["atom_grid"] == [75, 302]
        assert numerics["grid"]["num_points"] > 0
        assert output.spec["aux_basis"] is None
        assert numerics["auxbasis"] == "def2-svp-jkfit"

    def test_def2_basis_brings_its_core_potential(self, run_pyscf_driver):
        """Gaussian and ORCA apply def2's 28-electron potential to iodine.

        Without it PySCF would run iodine all-electron in a basis built for
        the potential, a different and meaningless calculation.
        """
        hydrogen_iodide = Molecule(
            symbols=["I", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.61]]),
            charge=0,
            multiplicity=1,
        )
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            hydrogen_iodide,
            "sp",
            **{**HF_MINIMAL, "basis": "def2-svp"},
        )
        assert output.ecp_core_electrons == {"I": 28}
        assert output.spec["num_electrons"] == 54 - 28
        mol = gto.M(
            atom="I 0 0 0; H 0 0 1.61",
            basis="def2-svp",
            ecp={"I": "def2-svp"},
            verbose=0,
        )
        assert output.energies[-1] == pytest.approx(
            hf.RHF(mol).kernel(), abs=1e-7
        )

    def test_hartree_fock_records_no_functional(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob, pyscf_h2_molecule, "sp", **HF_MINIMAL
        )
        assert output.reference == "RHF"
        assert output.spin == "restricted"
        assert output.xc_applied is None
        assert output.spin_square == 0.0
        assert output.numerics_applied["grid"] is None
        assert output.numerics_applied["auxbasis"] is None
        assert output.ecp_core_electrons is None

    def test_open_shell_reference_records_its_spin_contamination(
        self, run_pyscf_driver, pyscf_oh_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_oh_molecule,
            "sp",
            multiplicity=2,
            **HF_MINIMAL,
        )
        assert output.reference == "UHF"
        assert output.spin == "unrestricted"
        # A UHF determinant has <S^2> >= S(S+1); a doublet is only slightly
        # contaminated by quartet character.
        assert 0.0 <= output.spin_contamination < 0.05

    def test_configuration_redefining_b3lyp_is_refused(
        self, run_pyscf_driver, pyscf_h2_molecule, monkeypatch
    ):
        monkeypatch.setattr(__config__, "B3LYP_WITH_VWN5", True, raising=False)
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional="b3lyp",
            **DFT_MINIMAL,
        )
        assert returncode == 1
        assert output.failure["stage"] == "initialization"
        assert "B3LYP_WITH_VWN5" in output.failure["message"]
        assert "energies" not in output.results

    @pytest.mark.parametrize(
        ("functional", "hint"),
        [
            # B3LYP minus VWN_RPA: a functional nobody meant.
            ("b3lyp-vwn3", None),
            # PySCF knows only the underscore spelling.
            ("lc-wpbe", "lc_wpbe"),
        ],
    )
    def test_dash_read_as_subtraction_is_refused(
        self, run_pyscf_driver, pyscf_h2_molecule, functional, hint
    ):
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            functional=functional,
            **DFT_MINIMAL,
        )
        assert returncode == 1
        assert output.failure["stage"] == "initialization"
        assert "subtraction" in output.failure["message"]
        if hint is not None:
            assert repr(hint) in output.failure["message"]
        assert "energies" not in output.results

    @pytest.mark.parametrize(
        ("log", "basis", "reference", "tolerance"),
        [
            # Gaussian 16 B.01 B3LYP/def2-QZVP; one atom, so GD3BJ adds 0.
            ("he.log", "def2-qzvp", -2.91512971456, 1e-6),
            # Gaussian 16 B.01 B3LYP/def2-SVP: SCF Done less the GD3BJ term,
            # the nuclear repulsion after minus before dispersion.
            (
                "co2.log",
                "def2-svp",
                -188.444679593 - (58.2365923328 - 58.2382352385),
                1e-5,
            ),
        ],
    )
    def test_b3lyp_reproduces_gaussian_b3lyp(
        self,
        run_pyscf_driver,
        gaussian_outputs_test_directory,
        log,
        basis,
        reference,
        tolerance,
    ):
        """Far inside the 7e-3 to 8e-2 Eh that b3lyp5 lies away here."""
        molecule = Molecule.from_filepath(
            os.path.join(gaussian_outputs_test_directory, log), index="-1"
        )
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            molecule,
            "sp",
            functional="b3lyp",
            basis=basis,
            density_fit=False,
            defgrid="defgrid3",
            scf_tol=1e-10,
        )
        assert output.energies[-1] == pytest.approx(reference, abs=tolerance)

    def test_b3lyp5_reproduces_orca_b3lyp(
        self, run_pyscf_driver, orca_outputs_directory
    ):
        """ORCA 6.1.0 B3LYP (VWN5)/def2-SVP on KOH, SCF energy without D3.

        ORCA ran with RIJCOSX and this run uses exact exchange, so they agree
        to about 1e-4 Eh rather than 1e-6; b3lyp misses by 0.1 Eh.
        """
        molecule = Molecule.from_filepath(
            os.path.join(orca_outputs_directory, "KOH.out"), index="-1"
        )
        reference = -675.52250804211144
        energies = {}
        for functional in ("b3lyp5", "b3lyp"):
            _, output = run_pyscf_driver(
                PySCFSinglePointJob,
                molecule,
                "sp",
                label=functional,
                functional=functional,
                basis="def2-svp",
                density_fit=False,
                defgrid="defgrid2",
                scf_tol=1e-10,
            )
            energies[functional] = output.energies[-1]
        assert energies["b3lyp5"] == pytest.approx(reference, abs=1e-3)
        assert abs(energies["b3lyp"] - reference) > 5e-2

    def test_hartree_fock_is_unaffected_by_that_configuration(
        self, run_pyscf_driver, pyscf_h2_molecule, monkeypatch
    ):
        monkeypatch.setattr(__config__, "B3LYP_WITH_VWN5", True, raising=False)
        returncode, output = run_pyscf_driver(
            PySCFSinglePointJob, pyscf_h2_molecule, "sp", **HF_MINIMAL
        )
        assert returncode == 0
        assert output.normal_termination


class TestCorrelatedMethods:
    """MP2, CCSD and CCSD(T) single points on the Hartree-Fock reference."""

    def test_ccsd_of_two_electrons_is_full_ci(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        """With two electrons CCSD is exact and the (T) correction vanishes."""
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_h2_molecule,
            "sp",
            ab_initio="ccsd(t)",
            basis="cc-pvdz",
            correlation={"conv_tol": 1e-10},
            **CORRELATED,
        )
        mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="cc-pvdz", verbose=0)
        exact = fci.FCI(hf.RHF(mol).run(conv_tol=1e-11)).kernel()[0]
        assert output.normal_termination
        assert output.energies[-1] == pytest.approx(exact, abs=1e-8)
        assert output.triples_correction == pytest.approx(0.0, abs=1e-10)
        assert output.status["stages"]["correlation"]["frozen_orbitals"] == 0

    @pytest.mark.parametrize(
        ("orca_output", "frozen_core", "frozen_orbitals"),
        [
            ("water_mp2_frozen_core.out", True, 1),
            ("water_mp2_all_electron.out", False, 0),
        ],
    )
    def test_mp2_reproduces_orca(
        self,
        run_pyscf_driver,
        orca_outputs_directory,
        orca_output,
        frozen_core,
        frozen_orbitals,
    ):
        """ORCA 6.1.1 MP2/def2-SVP with exact integrals (NoRI)."""
        path = os.path.join(orca_outputs_directory, orca_output)
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            Molecule.from_filepath(path, index="-1"),
            "sp",
            ab_initio="mp2",
            basis="def2-svp",
            correlation={"frozen_core": frozen_core},
            **CORRELATED,
        )
        stage = output.status["stages"]["correlation"]
        assert output.energies[-1] == pytest.approx(
            ORCAOutput(path).energies[-1], abs=1e-6
        )
        assert stage["frozen_orbitals"] == frozen_orbitals
        assert stage["solver"] == "RMP2"

    @pytest.mark.parametrize(
        ("atom", "frozen"),
        [
            # ORCA 6.1's default frozen core (manual Table 2.69), in
            # orbitals, beyond any effective core potential.
            ("C 0 0 0; O 0 0 1.13", 2),
            ("H 0 0 0; Cl 0 0 1.27", 5),
            ("K 0 0 0; Cl 0 0 2.67", 10),
            ("Zn 0 0 0; Cl 0 0 2.05; Cl 0 0 -2.05", 15),
            ("H 0 0 0; Br 0 0 1.41", 9),
            # def2-SVP replaces iodine's 28 innermost electrons.
            ("H 0 0 0; I 0 0 1.61", 4),
        ],
    )
    def test_default_frozen_core_is_orcas(self, atom, frozen):
        """Fails loudly if a newer PySCF redefines its chemical core.

        Gaussian's FC differs for Li and Be, Sc-Zn, and Ga-Kr with Pople
        basis sets, so ORCA's table is the reference.
        """
        mol = gto.M(atom=atom, basis="def2-svp", ecp="def2-svp", verbose=0)
        assert mp.MP2(hf.RHF(mol)).set_frozen().frozen == frozen

    def test_recorded_frozen_core_is_the_one_applied(self, run_pyscf_driver):
        hydrogen_bromide = Molecule(
            symbols=["H", "Br"],
            positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.41]]),
            charge=0,
            multiplicity=1,
        )
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            hydrogen_bromide,
            "sp",
            ab_initio="mp2",
            basis="def2-svp",
            **CORRELATED,
        )
        assert output.status["stages"]["correlation"]["frozen_orbitals"] == 9

    def test_mp2_reproduces_gaussian(
        self, run_pyscf_driver, gaussian_outputs_test_directory
    ):
        """Gaussian 16 B.01 MP2/aug-cc-pVTZ with its default frozen core."""
        path = os.path.join(gaussian_outputs_test_directory, "water_mp2.log")
        reference = Gaussian16Output(path)
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            Molecule.from_filepath(path, index="-1"),
            "sp",
            ab_initio="mp2",
            basis="aug-cc-pvtz",
            **CORRELATED,
        )
        assert output.scf_energies[-1] == pytest.approx(
            reference.scf_energies[-1], abs=1e-6
        )
        assert output.energies[-1] == pytest.approx(
            reference.energies[-1], abs=1e-6
        )

    @pytest.mark.parametrize("method", ["ccsd", "ccsd(t)"])
    def test_coupled_cluster_matches_pyscf_run_directly(
        self, run_pyscf_driver, pyscf_water_molecule, method
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            ab_initio=method,
            basis="cc-pvdz",
            correlation={"conv_tol": 1e-10},
            **CORRELATED,
        )
        mol = gto.M(
            atom=list(
                zip(
                    pyscf_water_molecule.chemical_symbols,
                    pyscf_water_molecule.positions,
                )
            ),
            basis="cc-pvdz",
            verbose=0,
        )
        solver = cc.CCSD(hf.RHF(mol).run(conv_tol=1e-11))
        solver.conv_tol = 1e-10
        solver.set_frozen()
        solver.kernel()
        expected = solver.e_tot
        if method == "ccsd(t)":
            expected += solver.ccsd_t()
        assert output.energies[-1] == pytest.approx(expected, abs=1e-8)
        assert output.correlation_energy == pytest.approx(
            expected - output.scf_energies[-1], abs=1e-8
        )
        # The SCF density would describe the reference, not this result.
        assert "dipole_moment" in output.status["skipped_properties"]

    def test_open_shell_uses_the_unrestricted_reference(
        self, run_pyscf_driver, pyscf_oh_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_oh_molecule,
            "sp",
            ab_initio="ccsd(t)",
            basis="cc-pvdz",
            multiplicity=2,
            **CORRELATED,
        )
        assert output.normal_termination
        assert output.reference == "UHF"
        assert output.status["stages"]["correlation"]["solver"] == "UCCSD"
        assert output.triples_correction < 0

    def test_unconverged_coupled_cluster_reports_no_energy(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            ab_initio="ccsd(t)",
            basis="cc-pvdz",
            correlation={"max_cycle": 1, "conv_tol": 1e-12},
            **CORRELATED,
        )
        assert output.normal_termination is False
        assert output.energies == []
        assert output.correlation_energy is None
        assert output.triples_correction is None
        assert len(output.scf_energies) == 1
        assert "correlation" in output.failure_reason

    @pytest.mark.parametrize("method", ["mp2", "ccsd"])
    def test_forces_are_the_correlated_gradient(
        self, run_pyscf_driver, pyscf_distorted_water_molecule, method
    ):
        """Analytic forces against a central difference of the energy."""
        molecule = pyscf_distorted_water_molecule
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            molecule,
            "sp",
            ab_initio=method,
            basis="6-31g",
            forces=True,
            correlation={"conv_tol": 1e-10} if method == "ccsd" else None,
            **CORRELATED,
        )
        bohr = 0.52917721092
        step = 1e-3  # Bohr

        def energy(positions):
            mol = gto.M(
                atom=list(zip(molecule.chemical_symbols, positions)),
                basis="6-31g",
                verbose=0,
            )
            reference = hf.RHF(mol).run(conv_tol=1e-12)
            solver = (
                mp.MP2(reference) if method == "mp2" else cc.CCSD(reference)
            )
            if method == "ccsd":
                solver.conv_tol = 1e-11
            solver.set_frozen()
            return solver.kernel()[0] + reference.e_tot

        numerical = np.zeros_like(molecule.positions)
        for atom in range(len(molecule.positions)):
            for axis in range(3):
                shifted = [molecule.positions.copy() for _ in range(2)]
                shifted[0][atom, axis] += step * bohr
                shifted[1][atom, axis] -= step * bohr
                gradient = (energy(shifted[0]) - energy(shifted[1])) / (
                    2 * step
                )
                numerical[atom, axis] = -gradient
        assert np.abs(output.forces - numerical).max() < 1e-5

    @pytest.mark.slow
    def test_mp2_optimum_reproduces_gaussian(
        self, run_pyscf_driver, gaussian_outputs_test_directory
    ):
        """Gaussian 16 B.01 MP2(FC)/aug-cc-pVTZ optimum from a displaced start."""
        path = os.path.join(gaussian_outputs_test_directory, "water_mp2.log")
        optimum = Molecule.from_filepath(path, index="-1")
        start = optimum.copy()
        start.positions = optimum.positions + np.array(
            [[0.0, 0.0, 0.0], [0.03, 0.0, 0.0], [0.0, -0.02, 0.0]]
        )
        _, output = run_pyscf_driver(
            PySCFOptJob,
            start,
            "opt",
            ab_initio="mp2",
            basis="aug-cc-pvtz",
            freq=False,
            **CORRELATED,
        )

        assert output.normal_termination
        assert output.status["stages"]["opt"]["optimizer_converged"]
        final = output.optimized_structure
        bond, angle = final.get_distance(1, 2), final.get_angle(2, 1, 3)
        reference_bond = optimum.get_distance(1, 2)
        reference_angle = optimum.get_angle(2, 1, 3)
        assert abs(bond - reference_bond) < 5e-4
        assert abs(angle - reference_angle) < 0.1
        assert output.energies[-1] == pytest.approx(
            Gaussian16Output(path).energies[-1], abs=1e-6
        )

    def test_record_carries_the_correlated_total(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            ab_initio="mp2",
            basis="cc-pvdz",
            **CORRELATED,
        )
        record = PySCFAssembler(filename=output.filename).assemble()
        assert record.meta["method"] == "mp2"
        assert record.results["total_energy"] == pytest.approx(
            output.energies[-1]
        )
        assert record.results["total_energy"] < output.scf_energies[-1]


class TestSCFStability:
    """Stability of the converged reference, recorded and never followed."""

    def test_equilibrium_water_is_stable_in_every_direction(
        self, run_pyscf_driver, pyscf_water_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_water_molecule,
            "sp",
            scf_stability=True,
            **{**HF_MINIMAL, "basis": "cc-pvdz"},
        )
        record = output.scf_stability
        assert record["external"]["question"] == "RHF/RKS -> UHF/UKS"
        for question in ("internal", "real_to_complex", "external"):
            assert record[question]["stable"] is True
            assert record[question]["lowest_eigenvalue"] > 0

    def test_stretched_h2_is_unstable_towards_unrestricted(
        self, run_pyscf_driver
    ):
        """Past the Coulson-Fischer point RHF is no longer the UHF minimum.

        The instability is reported and the requested RHF solution kept.
        """
        geometry = (["H", "H"], [[0, 0, 0], [0, 0, 2.5]], 1)
        analysed = run_pyscf_ccpvdz_sp(
            run_pyscf_driver,
            *geometry,
            label="analysed",
            ab_initio="hf",
            scf_stability=True,
        )
        plain = run_pyscf_ccpvdz_sp(
            run_pyscf_driver, *geometry, label="plain", ab_initio="hf"
        )
        record = analysed.scf_stability
        assert record["internal"]["stable"] is True
        assert record["external"]["stable"] is False
        assert record["external"]["lowest_eigenvalue"] < -1e-5
        assert analysed.energies == plain.energies
        assert analysed.normal_termination

    def test_real_to_complex_instability_is_reported(self, run_pyscf_driver):
        """PySCF computes this answer but returns only the RKS -> UKS one."""
        output = run_pyscf_ccpvdz_sp(
            run_pyscf_driver,
            ["O", "O"],
            [[0, 0, 0], [0, 0, 1.21]],
            1,
            functional="b3lyp",
            scf_stability=True,
        )
        assert output.scf_stability["real_to_complex"]["stable"] is False

    def test_unrestricted_reference_asks_about_ghf(
        self, run_pyscf_driver, pyscf_oh_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob,
            pyscf_oh_molecule,
            "sp",
            multiplicity=2,
            scf_stability=True,
            **HF_MINIMAL,
        )
        assert output.scf_stability["reference"] == "UHF"
        assert (
            output.scf_stability["external"]["question"]
            == "UHF/UKS -> GHF/GKS"
        )

    def test_not_requested_is_not_reported_as_stable(
        self, run_pyscf_driver, pyscf_h2_molecule
    ):
        _, output = run_pyscf_driver(
            PySCFSinglePointJob, pyscf_h2_molecule, "sp", **HF_MINIMAL
        )
        assert output.spec["scf_stability"] is False
        assert output.scf_stability is None


class TestTransitionState:
    """HCN <=> HNC isomerisation at RHF/6-31G*.

    Three atoms, so the saddle point and its Hessian take seconds, and the
    reaction coordinate -- the hydrogen swinging from carbon to nitrogen --
    is unambiguous.
    """

    BASIS = {"ab_initio": "hf", "basis": "6-31g*", "density_fit": False}

    @classmethod
    def _search(cls, run_pyscf_driver, hydrogen, label, freq=True):
        molecule = Molecule(
            symbols=["C", "N", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.18], hydrogen]),
            charge=0,
            multiplicity=1,
        )
        return run_pyscf_driver(
            PySCFTSJob,
            molecule,
            "ts",
            label=label,
            freq=freq,
            scf_tol=1e-10,
            **cls.BASIS,
        )[1]

    def test_bridging_guess_converges_to_the_saddle_point(
        self, run_pyscf_driver
    ):
        output = self._search(run_pyscf_driver, [1.109, 0.0, 0.404], "ts")
        assert output.normal_termination
        assert output.jobtype == "ts"
        assert output.saddle_order == 1
        imaginary = output.status["stages"]["ts"]["imaginary_frequencies"]
        assert imaginary == pytest.approx([-1226.7], abs=1.0)
        assert output.vibrational_frequencies[1:] == pytest.approx(
            [2238.4, 2857.9], abs=1.0
        )
        # Checked independently of the optimiser: PySCF's own gradient
        # vanishes at the reported structure, whose energy it reproduces.
        coordinates = output.optimized_structure.positions
        mol = gto.M(
            atom=list(zip(["C", "N", "H"], coordinates.tolist())),
            basis="6-31g*",
            verbose=0,
        )
        mf = hf.RHF(mol).run(conv_tol=1e-10)
        assert mf.e_tot == pytest.approx(output.energies[-1], abs=1e-8)
        assert np.abs(mf.nuc_grad_method().kernel()).max() < 5e-4
        assert output.energies[-1] == pytest.approx(-92.791419, abs=1e-5)
        # The reaction coordinate is not a vibration of the transition
        # state: the zero-point energy counts only the two real modes.
        real = output.vibrational_frequencies[1:]
        hartree_per_wavenumber = 1 / 219474.6313632
        assert output.zero_point_energy == pytest.approx(
            0.5 * sum(real) * hartree_per_wavenumber, rel=1e-4
        )

    def test_minimum_found_instead_is_not_reported_as_a_transition_state(
        self, run_pyscf_driver
    ):
        """From a guess near HCN the search falls into the HCN minimum."""
        output = self._search(run_pyscf_driver, [1.05, 0.0, 0.75], "minimum")
        assert not output.normal_termination
        assert output.saddle_order == 0
        assert output.status["stages"]["ts"]["imaginary_frequencies"] == []
        assert "0 imaginary modes instead of one" in output.failure_reason
        # The structure and its frequencies are kept for inspection.
        assert min(output.vibrational_frequencies) > 0

    @pytest.mark.parametrize(
        ("frequencies", "order"),
        [
            ([-1226.7, 2238.4, 2857.9], 1),
            # No threshold: a small negative mode counts, as it does in
            # CHEMSMART's thermochemistry, so a second-order saddle is never
            # reported as a transition state.
            ([-1226.7, -3.0, 2857.9], 2),
            ([531.5, 2308.4, 4092.2], 0),
        ],
    )
    def test_every_negative_frequency_is_an_imaginary_mode(
        self, frequencies, order
    ):
        status = {"converged": True}
        driver._check_saddle_order(status, frequencies)
        assert status["saddle_order"] == order
        assert status["converged"] is (order == 1)

    def test_saddle_order_is_unknown_without_a_hessian(self, run_pyscf_driver):
        output = self._search(
            run_pyscf_driver, [1.109, 0.0, 0.404], "no_freq", freq=False
        )
        assert output.normal_termination
        assert output.spec["stages"] == ["scf", "ts"]
        assert output.saddle_order is None
        assert output.vibrational_frequencies is None


class TestIntrinsicReactionCoordinate:
    """Both branches from the HCN <=> HNC saddle point at RHF/6-31G*."""

    BASIS = {"ab_initio": "hf", "basis": "6-31g*", "density_fit": False}
    SYMBOLS = ["C", "N", "H"]
    #: The saddle point found by TestTransitionState, in Angstrom.
    SADDLE_POINT = [[0.0065, 0.0, 0.0446], [-0.0248, 0.0, 1.2132]]
    SADDLE_POINT += [[1.1267, 0.0, 0.3259]]

    @classmethod
    def _irc(cls, run_pyscf_driver, positions, label, **settings):
        molecule = Molecule(
            symbols=cls.SYMBOLS,
            positions=np.array(positions, dtype=float),
            charge=0,
            multiplicity=1,
        )
        return run_pyscf_driver(
            PySCFIRCJob,
            molecule,
            "irc",
            label=label,
            scf_tol=1e-10,
            **cls.BASIS,
            **settings,
        )[1]

    @staticmethod
    def _bonds(positions):
        """Return r(CN) and the hydrogen's distances to C and N, in A."""
        c, n, h = np.asarray(positions, dtype=float)
        return np.array(
            [
                np.linalg.norm(n - c),
                np.linalg.norm(h - c),
                np.linalg.norm(h - n),
            ]
        )

    @classmethod
    def _minimum(cls, positions):
        """Optimise a minimum with PySCF and geomeTRIC, outside CHEMSMART."""
        from pyscf.geomopt import geometric_solver

        mol = gto.M(
            atom=list(zip(cls.SYMBOLS, positions)),
            basis="6-31g*",
            verbose=0,
        )
        mf = hf.RHF(mol)
        mf.conv_tol = 1e-10
        converged, optimised = geometric_solver.kernel(mf, maxsteps=100)
        assert converged
        energy = hf.RHF(optimised).run(conv_tol=1e-10).e_tot
        return energy, cls._bonds(optimised.atom_coords(unit="Angstrom"))

    def test_branches_end_in_hydrogen_cyanide_and_isocyanide(
        self, run_pyscf_driver
    ):
        output = self._irc(run_pyscf_driver, self.SADDLE_POINT, "irc")
        assert output.normal_termination, output.failure_reason
        stage = output.status["stages"]["irc"]
        assert stage["saddle_order"] == 1
        assert stage["imaginary_frequencies"] == pytest.approx(
            [-1226.7], abs=1.0
        )
        assert all(b["converged"] for b in stage["branches"].values())

        path = output.irc_path
        start = output.irc_ts_index
        assert 0 < start < len(path) - 1
        assert np.array_equal(path[start].positions, output.positions)
        assert path[start].energy == output.energies[-1]
        energies = np.array([frame.energy for frame in path])
        # Downhill from the saddle point in both directions.
        assert np.all(np.diff(energies[: start + 1]) >= -1e-6)
        assert np.all(np.diff(energies[start:]) <= 1e-6)

        # The ends are the two minima, whichever branch reached which.
        hcn = self._minimum([[0, 0, 0], [0, 0, 1.15], [0, 0, -1.07]])
        hnc = self._minimum([[0, 0, 1.17], [0, 0, 0], [0, 0, -1.0]])
        ends = sorted(
            (path[0], path[-1]),
            key=lambda frame: self._bonds(frame.positions)[1],
        )
        for end, (energy, bonds) in zip(ends, (hcn, hnc)):
            assert end.energy == pytest.approx(energy, abs=1e-5)
            assert self._bonds(end.positions) == pytest.approx(bonds, abs=2e-3)

    def test_one_branch_is_the_same_branch_of_both(self, run_pyscf_driver):
        both = self._irc(run_pyscf_driver, self.SADDLE_POINT, "both")
        forward = self._irc(
            run_pyscf_driver,
            self.SADDLE_POINT,
            "forward",
            irc_direction="forward",
        )
        assert forward.normal_termination, forward.failure_reason
        assert forward.spec["irc_direction"] == "forward"
        assert list(forward.status["stages"]["irc"]["branches"]) == ["forward"]
        assert forward.irc_ts_index == 0
        assert np.abs(
            forward.irc_path[-1].positions - both.irc_path[-1].positions
        ).max() == pytest.approx(0.0, abs=1e-6)

    def test_branch_out_of_steps_is_not_reported_as_a_path(
        self, run_pyscf_driver
    ):
        output = self._irc(
            run_pyscf_driver, self.SADDLE_POINT, "short", opt_maxsteps=3
        )
        assert not output.normal_termination
        branches = output.status["stages"]["irc"]["branches"]
        assert [b["converged"] for b in branches.values()] == [False, False]
        assert "forward: no minimum within 3 steps" in output.failure_reason
        assert output.irc_path is None

    def test_minimum_is_not_followed_as_a_transition_state(
        self, run_pyscf_driver
    ):
        hydrogen_cyanide = [[0, 0, 0], [0, 0, 1.133], [0, 0, -1.059]]
        output = self._irc(run_pyscf_driver, hydrogen_cyanide, "minimum")
        assert not output.normal_termination
        assert output.status["stages"]["irc"]["saddle_order"] == 0
        assert "0 imaginary modes instead of one" in output.failure_reason
        assert output.irc_path is None
        # The start's Hessian is a check, not a frequency result.
        assert output.vibrational_frequencies is None
