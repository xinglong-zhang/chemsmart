"""CREST data-file and settings fixtures."""

import os

import pytest


# master CREST test directory
@pytest.fixture()
def crest_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "CRESTTests")


# crest yaml files
@pytest.fixture()
def crest_yaml_settings_directory(crest_test_directory):
    return os.path.join(crest_test_directory, "project_yaml")


@pytest.fixture()
def crest_yaml_settings_gas(crest_yaml_settings_directory):
    return os.path.join(crest_yaml_settings_directory, "gas.yaml")


@pytest.fixture()
def crest_yaml_settings_gas_project_name(crest_yaml_settings_directory):
    return os.path.join(crest_yaml_settings_directory, "gas")


@pytest.fixture()
def crest_yaml_settings_solv(crest_yaml_settings_directory):
    return os.path.join(crest_yaml_settings_directory, "solv.yaml")


@pytest.fixture()
def crest_yaml_settings_solv_project_name(crest_yaml_settings_directory):
    return os.path.join(crest_yaml_settings_directory, "solv")


@pytest.fixture()
def crest_outputs_directory(crest_test_directory):
    crest_outputs_directory = os.path.join(crest_test_directory, "outputs")
    return os.path.abspath(crest_outputs_directory)


@pytest.fixture()
def crest_styrene_outfolder(crest_outputs_directory):
    return os.path.join(crest_outputs_directory, "styrene_conformers")


@pytest.fixture()
def crest_octane_outfolder(crest_outputs_directory):
    return os.path.join(crest_outputs_directory, "octane_conformers")


@pytest.fixture()
def crest_ts1a_constrained_outfolder(crest_outputs_directory):
    return os.path.join(crest_outputs_directory, "TS1A_constrained_conformers")


@pytest.fixture()
def crest_hexane_failed_outfolder(crest_outputs_directory):
    return os.path.join(crest_outputs_directory, "hexane_failed_conformers")
