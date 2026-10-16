"""
Direct unit tests for :class:`AtomsChargeMultiplicity` in
``chemsmart.io.molecules.atoms``.
"""

import numpy as np
import pytest
from ase import Atoms
from ase.constraints import FixAtoms

from chemsmart.io.molecules.atoms import AtomsChargeMultiplicity


@pytest.fixture()
def simple_ase_atoms():
    return Atoms("Ar2", positions=[(0.0, 0.0, 0.0), (3.5, 0.0, 0.0)])


class TestChargeSetter:
    def test_accepts_int(self):
        obj = AtomsChargeMultiplicity(
            charge=0,
            multiplicity=1,
            frozen_atoms=None,
            energy=None,
            forces=None,
            symbols=["Ar"],
            positions=[(0.0, 0.0, 0.0)],
        )
        obj.charge = 1
        assert obj.charge == 1

    def test_accepts_integer_valued_float(self):
        obj = AtomsChargeMultiplicity(
            charge=0,
            multiplicity=1,
            frozen_atoms=None,
            energy=None,
            forces=None,
            symbols=["Ar"],
            positions=[(0.0, 0.0, 0.0)],
        )
        obj.charge = 2.0
        assert obj.charge == 2
        assert isinstance(obj.charge, int)

    def test_rejects_non_integer_float(self):
        obj = AtomsChargeMultiplicity(
            charge=0,
            multiplicity=1,
            frozen_atoms=None,
            energy=None,
            forces=None,
            symbols=["Ar"],
            positions=[(0.0, 0.0, 0.0)],
        )
        with pytest.raises(TypeError, match="Charge must be an integer"):
            obj.charge = 1.5

    def test_rejects_string(self):
        obj = AtomsChargeMultiplicity(
            charge=0,
            multiplicity=1,
            frozen_atoms=None,
            energy=None,
            forces=None,
            symbols=["Ar"],
            positions=[(0.0, 0.0, 0.0)],
        )
        with pytest.raises(TypeError, match="Charge must be an integer"):
            obj.charge = "1"


class TestMultiplicitySetter:
    def test_accepts_int(self):
        obj = AtomsChargeMultiplicity(
            charge=0,
            multiplicity=1,
            frozen_atoms=None,
            energy=None,
            forces=None,
            symbols=["Ar"],
            positions=[(0.0, 0.0, 0.0)],
        )
        obj.multiplicity = 3
        assert obj.multiplicity == 3

    def test_rejects_non_integer_float(self):
        obj = AtomsChargeMultiplicity(
            charge=0,
            multiplicity=1,
            frozen_atoms=None,
            energy=None,
            forces=None,
            symbols=["Ar"],
            positions=[(0.0, 0.0, 0.0)],
        )
        with pytest.raises(
            TypeError, match="Spin multiplicity must be an integer"
        ):
            obj.multiplicity = 2.5


class TestFromAtomsNoConstraints:
    def test_no_constraints_gives_none_frozen_atoms(self, simple_ase_atoms):
        result = AtomsChargeMultiplicity.from_atoms(
            simple_ase_atoms, charge=0, multiplicity=1
        )
        assert result.frozen_atoms is None

    def test_no_calculator_leaves_energy_forces_velocities_none(
        self, simple_ase_atoms
    ):
        result = AtomsChargeMultiplicity.from_atoms(
            simple_ase_atoms, charge=0, multiplicity=1
        )
        assert result.energy is None
        assert result.forces is None


class TestFromAtomsSingleFixAtomsConstraint:
    def test_single_fixatoms_constraint_not_wrapped_in_list(
        self, simple_ase_atoms
    ):
        """The `elif isinstance(atoms.constraints, FixAtoms):` branch in
        from_atoms is dead code through any public ASE API: ASE's own
        `Atoms.constraints` setter (aliased to set_constraint) always
        wraps a bare constraint in a list
        (`self._constraints = [constraint]`), so `atoms.constraints`
        can never actually be a bare FixAtoms instance -- see
        BUGS_FOUND.md. The only way to reach this branch at all is by
        writing the private `_constraints` attribute directly,
        bypassing the public property entirely, as done here.
        """
        simple_ase_atoms._constraints = FixAtoms(indices=[0])
        result = AtomsChargeMultiplicity.from_atoms(
            simple_ase_atoms, charge=0, multiplicity=1
        )
        assert result.frozen_atoms == [-1, 0]


class TestFromAtomsCalculatorError:
    def test_calculator_error_wrapped_as_runtime_error(self, simple_ase_atoms):
        from ase.calculators.calculator import Calculator, CalculatorError

        class _BrokenCalculator(Calculator):
            implemented_properties = ["energy", "forces"]

            def calculate(self, *args, **kwargs):
                raise CalculatorError("boom")

        simple_ase_atoms.calc = _BrokenCalculator()
        with pytest.raises(RuntimeError, match="Failed to obtain energy"):
            AtomsChargeMultiplicity.from_atoms(
                simple_ase_atoms, charge=0, multiplicity=1
            )


class TestToMolecule:
    def test_uses_own_charge_and_multiplicity_when_not_given(
        self, simple_ase_atoms
    ):
        obj = AtomsChargeMultiplicity.from_atoms(
            simple_ase_atoms, charge=1, multiplicity=2
        )
        mol = obj.to_molecule()
        assert mol.charge == 1
        assert mol.multiplicity == 2
        assert np.all(mol.symbols == ["Ar", "Ar"])
