import numpy as np
from ase import units

from chemsmart.utils.constants import amu_to_kg, bohr_to_meter


class TestUnits:

    def test_units(self):
        assert bohr_to_meter == 1 * units.Bohr / units.m
        assert np.isclose(bohr_to_meter, 0.52917721067e-10, atol=1e-10)
        assert amu_to_kg == 1 * units._amu
        assert np.isclose(amu_to_kg, 1.66053906660e-27, atol=1e-27)
