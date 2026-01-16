"""
Thermochemistry jobs module for quantum chemistry calculations.

This module provides comprehensive thermochemical analysis tools for
calculating thermal properties, entropies, enthalpies, and Gibbs free
energies from quantum chemistry output files.
"""

from .boltzmann import BoltzmannAverageThermochemistryJob
from .job import ThermochemistryJob
from .runner import ThermochemistryJobRunner

# Dynamically collect all thermochemistry job subclasses
jobs = ThermochemistryJob.subclasses()

__all__ = [
    "BoltzmannAverageThermochemistryJob",
    "ThermochemistryJob",
    "ThermochemistryJobRunner",
    "jobs",
]
