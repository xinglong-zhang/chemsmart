import numpy as np
import pytest
from ase import units

from chemsmart.utils.constants import (
    amu_to_kg,
    bohr_to_meter,
    energy_conversion,
    hartree_to_joules,
    joule_per_mol_to_eV,
)

# Define conversion factors for reference (in J/mol)
HARTREE_TO_J_MOL = hartree_to_joules * units._Nav  # ~2625.499639 kJ/mol
EV_TO_J_MOL = 1 / joule_per_mol_to_eV  # ~96.485339 kJ/mol
KCAL_MOL_TO_J_MOL = 4184.0  # kcal/mol to J/mol
KJ_MOL_TO_J_MOL = 1000.0  # kJ/mol to J/mol

# Tolerance for floating-point comparisons
TOLERANCE = 1e-6


class TestUnits:

    def test_units(self):
        assert bohr_to_meter == 1 * units.Bohr / units.m
        assert np.isclose(bohr_to_meter, 0.52917721067e-10, atol=1e-10)
        assert amu_to_kg == 1 * units._amu
        assert np.isclose(amu_to_kg, 1.66053906660e-27, atol=1e-27)


class TestEnergyConversion:
    def test_energy_conversion_hartree_to_kcal_mol(self):
        """Test conversion from Hartree to kcal/mol."""
        expected = HARTREE_TO_J_MOL / KCAL_MOL_TO_J_MOL  # ~627.509468
        result = energy_conversion("hartree", "kcal/mol")
        assert (
            abs(result - expected) < TOLERANCE
        ), f"Expected {expected}, got {result}"

    def test_energy_conversion_kcal_mol_to_hartree(self):
        """Test conversion from kcal/mol to Hartree."""
        value = 627.509468
        expected = value * KCAL_MOL_TO_J_MOL / HARTREE_TO_J_MOL  # ~1.0
        result = energy_conversion("kcal/mol", "hartree", value)
        assert (
            abs(result - expected) < TOLERANCE
        ), f"Expected {expected}, got {result}"

    def test_energy_conversion_ev_to_kj_mol(self):
        """Test conversion from eV to kJ/mol."""
        value = 1.0
        expected = EV_TO_J_MOL / KJ_MOL_TO_J_MOL  # ~96.485339
        result = energy_conversion("eV", "kJ/mol", value)
        assert (
            abs(result - expected) < TOLERANCE
        ), f"Expected {expected}, got {result}"

    def test_energy_conversion_j_mol_to_ev(self):
        """Test conversion from J/mol to eV."""
        value = 1000.0
        expected = value / EV_TO_J_MOL  # ~0.010364
        result = energy_conversion("J/mol", "eV", value)
        assert (
            abs(result - expected) < TOLERANCE
        ), f"Expected {expected}, got {result}"

    def test_energy_conversion_same_unit(self):
        """Test conversion with same unit (should return input value)."""
        value = 42.0
        for unit in ["hartree", "eV", "kcal/mol", "kJ/mol", "J/mol"]:
            result = energy_conversion(unit, unit, value)
            assert (
                abs(result - value) < TOLERANCE
            ), f"Expected {value}, got {result} for {unit}"

    def test_energy_conversion_zero_value(self):
        """Test conversion with zero input value."""
        result = energy_conversion("hartree", "kcal/mol", 0.0)
        assert abs(result - 0.0) < TOLERANCE, f"Expected 0.0, got {result}"

    def test_energy_conversion_negative_value(self):
        """Test conversion with negative input value."""
        value = -1.0
        expected = -HARTREE_TO_J_MOL / KCAL_MOL_TO_J_MOL  # ~-627.509468
        result = energy_conversion("hartree", "kcal/mol", value)
        assert (
            abs(result - expected) < TOLERANCE
        ), f"Expected {expected}, got {result}"

    def test_energy_conversion_invalid_from_unit(self):
        """Test error handling for invalid from_unit."""
        with pytest.raises(
            ValueError, match="Unsupported from_unit: invalid. Choose from"
        ):
            energy_conversion("invalid", "kcal/mol", 1.0)

    def test_energy_conversion_invalid_to_unit(self):
        """Test error handling for invalid to_unit."""
        with pytest.raises(
            ValueError, match="Unsupported to_unit: invalid. Choose from"
        ):
            energy_conversion("hartree", "invalid", 1.0)

    def test_energy_conversion_case_insensitivity(self):
        """Test case insensitivity for unit names."""
        value = 1.0
        expected = HARTREE_TO_J_MOL / KCAL_MOL_TO_J_MOL  # ~627.509468
        result = energy_conversion("HARTREE", "KCAL/MOL", value)
        assert (
            abs(result - expected) < TOLERANCE
        ), f"Expected {expected}, got {result}"

    def test_energy_conversion_logging(self, capture_log):
        """Test that logging captures conversion details."""
        value = 1.0
        energy_conversion("hartree", "kcal/mol", value)
        log_output = capture_log.text
        assert (
            "Converted 1.0 hartree to" in log_output
        ), "Logging message not found"
        assert "kcal/mol" in log_output, "Target unit not logged"

    def test_energy_conversion_all_pairs(self):
        """Test all possible unit conversion pairs."""
        value = 1.0
        units = ["hartree", "eV", "kcal/mol", "kJ/mol", "J/mol"]
        to_j_per_mol = {
            "hartree": HARTREE_TO_J_MOL,
            "eV": EV_TO_J_MOL,
            "kcal/mol": KCAL_MOL_TO_J_MOL,
            "kJ/mol": KJ_MOL_TO_J_MOL,
            "J/mol": 1.0,
        }
        for from_unit in units:
            for to_unit in units:
                expected = (
                    value * to_j_per_mol[from_unit] / to_j_per_mol[to_unit]
                )
                result = energy_conversion(from_unit, to_unit, value)
                assert (
                    abs(result - expected) < TOLERANCE
                ), f"Conversion from {from_unit} to {to_unit}: expected {expected}, got {result}"
