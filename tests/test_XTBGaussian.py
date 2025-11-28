"""
Tests for XTB calculator and Gaussian XTB job interface.

This module tests the XTBCalculator class and the Gaussian job
classes that use xtb as an external calculator for energy and
gradient calculations.
"""

import pytest

from chemsmart.calculators.xtb import XTBCalculator
from chemsmart.jobs.gaussian.xtb import (
    GaussianXTBJob,
    GaussianXTBJobSettings,
    GaussianXTBTSJob,
)


class TestXTBCalculator:
    """Tests for XTBCalculator configuration class."""

    def test_default_calculator(self):
        """Test default calculator initialization."""
        calc = XTBCalculator.default()
        assert calc.method == "gfn2"
        assert calc.num_threads == 1
        assert calc.solvent is None
        assert calc.solvation is None

    def test_calculator_with_threads(self):
        """Test calculator with parallel threads."""
        calc = XTBCalculator(num_threads=4)
        assert calc.num_threads == 4
        assert "-P 4" in calc.command_string

    def test_calculator_with_solvent(self):
        """Test calculator with implicit solvation."""
        calc = XTBCalculator(solvent="water", solvation="alpb")
        assert calc.solvent == "water"
        assert calc.solvation == "alpb"
        assert "--alpb water" in calc.command_string

    def test_calculator_with_gbsa(self):
        """Test calculator with GBSA solvation model."""
        calc = XTBCalculator(solvent="methanol", solvation="gbsa")
        assert "--gbsa methanol" in calc.command_string

    def test_calculator_with_electronic_temperature(self):
        """Test calculator with electronic temperature."""
        calc = XTBCalculator(electronic_temperature=300.0)
        assert "--etemp 300.0" in calc.command_string

    def test_calculator_with_method(self):
        """Test calculator with non-default method."""
        calc = XTBCalculator(method="gfn1")
        assert calc.method == "gfn1"
        assert "--gfn1" in calc.command_string

    def test_calculator_with_uhf(self):
        """Test calculator with unpaired electrons (spin state)."""
        calc = XTBCalculator(uhf=2)
        assert "--uhf 2" in calc.command_string

    def test_calculator_command_args(self):
        """Test command arguments list generation."""
        calc = XTBCalculator(
            num_threads=4,
            solvent="water",
            electronic_temperature=300.0,
        )
        args = calc.command_args
        assert "-P" in args
        assert "4" in args
        assert "--alpb" in args
        assert "water" in args
        assert "--etemp" in args
        assert "300.0" in args

    def test_external_command(self):
        """Test external command generation for Gaussian."""
        calc = XTBCalculator(num_threads=4, solvent="methanol")
        cmd = calc.external_command
        assert "xtb-gaussian" in cmd
        assert "-P 4" in cmd
        assert "--alpb methanol" in cmd

    def test_calculator_copy(self):
        """Test calculator deep copy."""
        calc = XTBCalculator(num_threads=8, solvent="water")
        calc_copy = calc.copy()
        assert calc_copy.num_threads == 8
        assert calc_copy.solvent == "water"
        # Modify copy and verify original is unchanged
        calc_copy.num_threads = 4
        assert calc.num_threads == 8

    def test_invalid_method_raises_error(self):
        """Test that invalid method raises ValueError."""
        with pytest.raises(ValueError, match="not supported"):
            XTBCalculator(method="invalid")

    def test_invalid_solvation_raises_error(self):
        """Test that invalid solvation model raises ValueError."""
        with pytest.raises(ValueError, match="not supported"):
            XTBCalculator(solvation="invalid")

    def test_calculator_to_dict(self):
        """Test conversion to dictionary."""
        calc = XTBCalculator(num_threads=4, solvent="water")
        d = calc.to_dict()
        assert d["num_threads"] == 4
        assert d["solvent"] == "water"
        assert d["method"] == "gfn2"

    def test_calculator_from_dict(self):
        """Test creation from dictionary."""
        config = {"num_threads": 4, "solvent": "methanol", "solvation": "gbsa"}
        calc = XTBCalculator.from_dict(config)
        assert calc.num_threads == 4
        assert calc.solvent == "methanol"
        assert calc.solvation == "gbsa"

    def test_log_all_flag(self):
        """Test log-all flag for verbose output."""
        calc = XTBCalculator(log_all=True)
        args = calc.command_args
        assert "--log-all" in args
        # log-all should be first argument
        assert args[0] == "--log-all"


class TestGaussianXTBJobSettings:
    """Tests for GaussianXTBJobSettings class."""

    def test_default_settings(self):
        """Test default XTB job settings."""
        settings = GaussianXTBJobSettings.default()
        assert settings.calculator is not None
        assert settings.job_type == "opt"
        assert settings.freq is False

    def test_settings_with_calculator(self):
        """Test settings with custom calculator."""
        calc = XTBCalculator(num_threads=8, solvent="water")
        settings = GaussianXTBJobSettings(calculator=calc)
        assert settings.calculator.num_threads == 8
        assert settings.calculator.solvent == "water"

    def test_settings_for_ts_optimization(self):
        """Test factory method for TS optimization."""
        settings = GaussianXTBJobSettings.for_ts_optimization(
            num_threads=4, solvent="methanol"
        )
        assert settings.job_type == "ts"
        assert settings.freq is True
        assert settings.calculator.num_threads == 4
        assert settings.calculator.solvent == "methanol"

    def test_settings_for_qst2(self):
        """Test factory method for QST2 calculation."""
        settings = GaussianXTBJobSettings.for_qst2(
            num_threads=4, solvent="water"
        )
        assert settings.job_type == "qst2"
        assert settings.freq is False
        assert settings.calculator.num_threads == 4

    def test_route_string_for_opt(self):
        """Test route string generation for optimization."""
        settings = GaussianXTBJobSettings(job_type="opt")
        route = settings.route_string
        assert "external=" in route
        assert "xtb-gaussian" in route
        assert "opt" in route

    def test_route_string_for_ts(self):
        """Test route string generation for TS optimization."""
        settings = GaussianXTBJobSettings(job_type="ts")
        route = settings.route_string
        assert "external=" in route
        assert "opt=(ts,calcfc,noeigentest)" in route

    def test_route_string_for_qst2(self):
        """Test route string generation for QST2."""
        settings = GaussianXTBJobSettings(job_type="qst2")
        route = settings.route_string
        assert "external=" in route
        assert "qst2" in route

    def test_route_string_with_freq(self):
        """Test route string includes frequency keyword."""
        settings = GaussianXTBJobSettings(job_type="opt", freq=True)
        route = settings.route_string
        assert "freq" in route

    def test_route_string_with_additional_options(self):
        """Test route string with additional optimization options."""
        settings = GaussianXTBJobSettings(
            job_type="ts",
            additional_opt_options_in_route="calcall,nomicro",
        )
        route = settings.route_string
        assert "calcall" in route
        assert "nomicro" in route

    def test_external_command_property(self):
        """Test external_command property."""
        calc = XTBCalculator(num_threads=4)
        settings = GaussianXTBJobSettings(calculator=calc)
        cmd = settings.external_command
        assert "-P 4" in cmd

    def test_settings_copy(self):
        """Test settings deep copy."""
        settings = GaussianXTBJobSettings(job_type="ts")
        settings_copy = settings.copy()
        settings_copy.job_type = "opt"
        assert settings.job_type == "ts"


class TestGaussianXTBJob:
    """Tests for GaussianXTBJob class."""

    @pytest.fixture
    def water_molecule(self):
        """Create a simple water molecule fixture."""
        from chemsmart.io.molecules.structure import Molecule

        symbols = ["O", "H", "H"]
        positions = [
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ]
        return Molecule(
            symbols=symbols,
            positions=positions,
            charge=0,
            multiplicity=1,
        )

    def test_job_initialization(self, water_molecule):
        """Test basic job initialization."""
        settings = GaussianXTBJobSettings.default()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianXTBJob(
            molecule=water_molecule,
            settings=settings,
            label="water_opt",
        )
        assert job.molecule is not None
        assert job.label == "water_opt"
        assert job.TYPE == "g16xtb"

    def test_job_calculator_property(self, water_molecule):
        """Test calculator property access."""
        calc = XTBCalculator(num_threads=4)
        settings = GaussianXTBJobSettings(calculator=calc, charge=0, multiplicity=1)
        job = GaussianXTBJob(
            molecule=water_molecule,
            settings=settings,
            label="test",
        )
        assert job.calculator.num_threads == 4


class TestGaussianXTBTSJob:
    """Tests for GaussianXTBTSJob class."""

    @pytest.fixture
    def ts_molecule(self):
        """Create a simple TS guess molecule fixture."""
        from chemsmart.io.molecules.structure import Molecule

        # Simple diatomic for testing
        symbols = ["H", "H"]
        positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        return Molecule(
            symbols=symbols,
            positions=positions,
            charge=0,
            multiplicity=1,
        )

    def test_ts_job_initialization(self, ts_molecule):
        """Test TS job initialization."""
        job = GaussianXTBTSJob(
            molecule=ts_molecule,
            label="ts_test",
            charge=0,
            multiplicity=1,
        )
        assert job.TYPE == "g16xtbts"
        assert job.settings.job_type == "ts"

    def test_ts_job_with_solvent(self, ts_molecule):
        """Test TS job with solvation."""
        settings = GaussianXTBJobSettings.for_ts_optimization(
            num_threads=4,
            solvent="water",
            charge=0,
            multiplicity=1,
        )
        job = GaussianXTBTSJob(
            molecule=ts_molecule,
            settings=settings,
            label="ts_solv",
        )
        assert job.calculator.solvent == "water"
