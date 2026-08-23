"""
Direct unit tests for chemsmart.jobs.gaussian.settings.GaussianpKaJobSettings
covering reference-acid validation/construction, conjugate-base/reference
charge-multiplicity fallback branches, and protonated-molecule handling
that lacked direct test coverage (tests/test_GaussianSettings.py covers
the main happy paths).
"""

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings


def _mol_with_h(charge=None, multiplicity=None, frozen_atoms=None):
    mol = Molecule(
        symbols=["C", "H", "H", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [0.63, -0.63, -0.63],
            [-0.63, 0.63, -0.63],
            [-0.63, -0.63, 0.63],
        ],
        frozen_atoms=frozen_atoms,
    )
    mol.charge = charge
    mol.multiplicity = multiplicity
    return mol


class TestBackwardsCompatThermodynamicCycle:
    def test_thermodynamic_cycle_kwarg_maps_to_scheme(self, caplog):
        with caplog.at_level("WARNING"):
            settings = GaussianpKaJobSettings(
                proton_index=2, thermodynamic_cycle="direct"
            )
        assert settings.scheme == "direct"
        assert "deprecated" in caplog.text


class TestTitleAndDelegatingMethods:
    def test_explicit_title_is_not_overridden(self):
        settings = GaussianpKaJobSettings(proton_index=2, title="My Title")
        assert settings.title == "My Title"

    def test_default_title_applied_when_falsy(self):
        settings = GaussianpKaJobSettings(proton_index=2, title=None)
        assert settings.title == "Gaussian pKa calculation job"

    def test_reference_pair_molecules_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
        )
        from unittest.mock import patch

        with patch.object(
            settings, "get_reference_molecule", return_value=_mol_with_h()
        ):
            acid, base = settings.reference_pair_molecules()
        assert acid is not None
        assert base.charge == -1

    def test_reference_pair_job_settings_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            functional="b3lyp",
            basis="sto-3g",
        )
        acid, base = settings.reference_pair_job_settings()
        assert acid.jobtype == "opt"
        assert base.jobtype == "opt"

    def test_reference_pair_sp_job_settings_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            functional="b3lyp",
            basis="sto-3g",
        )
        acid, base = settings.reference_pair_sp_job_settings()
        assert acid.jobtype == "sp"
        assert base.jobtype == "sp"

    def test_conjugate_pair_sp_job_settings_delegates(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, functional="b3lyp", basis="sto-3g"
        )
        mol = _mol_with_h(charge=0, multiplicity=1)
        prot, base = settings.conjugate_pair_sp_job_settings(mol)
        assert prot.jobtype == "sp"
        assert base.jobtype == "sp"


class TestValidateReferenceSettings:
    def test_no_reference_file_is_a_no_op(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        settings.validate_reference_settings()  # should not raise

    def test_missing_reference_charge_raises(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=1,
            reference_multiplicity=1,
        )
        with pytest.raises(ValueError, match="reference_charge"):
            settings.validate_reference_settings()

    def test_missing_reference_multiplicity_raises(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=1,
            reference_charge=0,
        )
        with pytest.raises(ValueError, match="reference_multiplicity"):
            settings.validate_reference_settings()


class TestGetReferenceMolecule:
    def test_raises_without_reference_file(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        with pytest.raises(ValueError, match="Reference file not provided"):
            settings.get_reference_molecule()


class TestCreateReferenceConjugateBaseMolecule:
    def _settings(self, **overrides):
        kwargs = dict(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
        )
        kwargs.update(overrides)
        return GaussianpKaJobSettings(**kwargs)

    def test_missing_reference_proton_index_raises(self):
        settings = self._settings(reference_proton_index=None)
        with pytest.raises(
            ValueError, match="reference_proton_index must be specified"
        ):
            settings._create_reference_conjugate_base_molecule(_mol_with_h())

    def test_out_of_range_reference_proton_index_raises(self):
        settings = self._settings(reference_proton_index=999)
        with pytest.raises(ValueError, match="out of range"):
            settings._create_reference_conjugate_base_molecule(_mol_with_h())

    def test_non_hydrogen_reference_proton_index_raises(self):
        settings = self._settings(reference_proton_index=1)  # index 1 = C
        with pytest.raises(ValueError, match="not hydrogen"):
            settings._create_reference_conjugate_base_molecule(_mol_with_h())

    def test_frozen_atoms_are_filtered(self):
        settings = self._settings(reference_proton_index=2)
        mol = _mol_with_h(frozen_atoms=[0, 0, 0, 0, 0])
        result = settings._create_reference_conjugate_base_molecule(mol)
        assert len(result.frozen_atoms) == 4

    def test_explicit_reference_conjugate_base_charge_and_multiplicity(self):
        settings = self._settings(
            reference_proton_index=2,
            reference_conjugate_base_charge=-5,
            reference_conjugate_base_multiplicity=3,
        )
        result = settings._create_reference_conjugate_base_molecule(
            _mol_with_h()
        )
        assert result.charge == -5
        assert result.multiplicity == 3

    def test_default_reference_conjugate_base_charge_and_multiplicity(self):
        settings = self._settings(reference_proton_index=2)
        result = settings._create_reference_conjugate_base_molecule(
            _mol_with_h()
        )
        assert result.charge == -1  # reference_charge - 1
        assert result.multiplicity == 1  # reference_multiplicity


class TestReferenceGasPhaseAndSolutionPhaseSettings:
    def _settings(self, **overrides):
        kwargs = dict(
            proton_index=2,
            reference_file="ref.xyz",
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            functional="b3lyp",
            basis="sto-3g",
        )
        kwargs.update(overrides)
        return GaussianpKaJobSettings(**kwargs)

    def test_gas_phase_default_conjugate_base_charge_and_mult(self):
        settings = self._settings()
        acid, base = settings._create_reference_gas_phase_job_settings()
        assert acid.charge == 0
        assert base.charge == -1
        assert base.multiplicity == 1

    def test_gas_phase_explicit_conjugate_base_charge_and_mult(self):
        settings = self._settings(
            reference_conjugate_base_charge=-7,
            reference_conjugate_base_multiplicity=2,
        )
        _, base = settings._create_reference_gas_phase_job_settings()
        assert base.charge == -7
        assert base.multiplicity == 2

    def test_solution_phase_default_conjugate_base_charge_and_mult(self):
        settings = self._settings()
        acid, base = settings._create_reference_solution_phase_sp_settings()
        assert acid.charge == 0
        assert base.charge == -1
        assert base.multiplicity == 1

    def test_solution_phase_explicit_conjugate_base_charge_and_mult(self):
        settings = self._settings(
            reference_conjugate_base_charge=-7,
            reference_conjugate_base_multiplicity=2,
        )
        _, base = settings._create_reference_solution_phase_sp_settings()
        assert base.charge == -7
        assert base.multiplicity == 2


class TestProtonatedMolecule:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, charge=2, multiplicity=3
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        result = settings.protonated_molecule(mol)
        assert result.charge == 2
        assert result.multiplicity == 3

    def test_falls_back_to_molecule_charge_and_multiplicity(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        settings.charge = None
        settings.multiplicity = None
        mol = _mol_with_h(charge=None, multiplicity=None)
        result = settings.protonated_molecule(mol)
        assert result.charge == 0
        assert result.multiplicity == 1


class TestCreateGasPhaseJobSettings:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=1,
            multiplicity=2,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_gas_phase_job_settings(mol)
        assert prot.charge == 1
        assert prot.multiplicity == 2

    def test_explicit_conjugate_base_charge_and_multiplicity(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=0,
            multiplicity=1,
            conjugate_base_charge=-9,
            conjugate_base_multiplicity=4,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        _, base = settings._create_gas_phase_job_settings(mol)
        assert base.charge == -9
        assert base.multiplicity == 4


class TestCreateMolecules:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, charge=3, multiplicity=2
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_molecules(mol)
        assert prot.charge == 3
        assert prot.multiplicity == 2

    def test_falls_back_to_defaults_when_molecule_charge_none(self):
        settings = GaussianpKaJobSettings(proton_index=2)
        settings.charge = None
        settings.multiplicity = None
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_molecules(mol)
        assert prot.charge == 0
        assert prot.multiplicity == 1


class TestCreateSolutionPhaseSpSettings:
    def test_explicit_settings_charge_and_multiplicity_used(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=1,
            multiplicity=2,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_solution_phase_sp_settings(mol)
        assert prot.charge == 1
        assert prot.multiplicity == 2

    def test_explicit_conjugate_base_charge_and_multiplicity(self):
        settings = GaussianpKaJobSettings(
            proton_index=2,
            functional="b3lyp",
            basis="sto-3g",
            charge=0,
            multiplicity=1,
            conjugate_base_charge=-9,
            conjugate_base_multiplicity=4,
        )
        mol = _mol_with_h(charge=None, multiplicity=None)
        _, base = settings._create_solution_phase_sp_settings(mol)
        assert base.charge == -9
        assert base.multiplicity == 4

    def test_falls_back_to_defaults_when_molecule_charge_none(self):
        settings = GaussianpKaJobSettings(
            proton_index=2, functional="b3lyp", basis="sto-3g"
        )
        settings.charge = None
        settings.multiplicity = None
        mol = _mol_with_h(charge=None, multiplicity=None)
        prot, _ = settings._create_solution_phase_sp_settings(mol)
        assert prot.charge == 0
        assert prot.multiplicity == 1
