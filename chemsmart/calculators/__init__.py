"""
Calculator interfaces for computational chemistry.

This module provides wrapper classes for interfacing various
quantum chemistry calculators with job management systems.
Enables use of external calculators like xtb with optimization
engines from packages like Gaussian and ORCA.
"""

from chemsmart.calculators.xtb import XTBCalculator

__all__ = ["XTBCalculator"]
