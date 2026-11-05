"""Thermochemistry fixtures."""

import shutil

import pytest


@pytest.fixture()
def gaussian_co2_he_output_directory(
    tmp_path, gaussian_co2_opt_outfile, gaussian_he_opt_outfile
):
    """``tmp_path`` holding copies of the Gaussian CO2 and He optimisation
    outputs (``co2.log``, ``he.log``), for directory-mode thermochemistry
    tests."""
    shutil.copy(gaussian_co2_opt_outfile, tmp_path / "co2.log")
    shutil.copy(gaussian_he_opt_outfile, tmp_path / "he.log")
    return tmp_path
