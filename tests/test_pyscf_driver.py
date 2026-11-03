"""Run the generated PySCF driver for real, on molecules small enough to be fast.

Skipped where PySCF is not installed, such as Windows, which it does not support.
Every number here comes from an actual calculation.
"""

from unittest.mock import patch

import numpy as np
import pytest

pytest.importorskip("pyscf")

from pyscf import gto  # noqa: E402
from pyscf import symm  # noqa: E402
from pyscf.scf import hf  # noqa: E402

from chemsmart.io.molecules.structure import Molecule  # noqa: E402
from chemsmart.io.pyscf import driver  # noqa: E402
from chemsmart.jobs.pyscf.hess import PySCFHessJob  # noqa: E402
from chemsmart.jobs.pyscf.opt import PySCFOptJob  # noqa: E402
from chemsmart.jobs.pyscf.settings import PySCFJobSettings  # noqa: E402
from chemsmart.jobs.pyscf.singlepoint import (  # noqa: E402
    PySCFSinglePointJob,
)

HF_MINIMAL = {"ab_initio": "hf", "basis": "sto-3g", "density_fit": False}


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
