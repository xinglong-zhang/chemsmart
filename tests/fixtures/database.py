"""Database fixtures."""

import os

import pytest


############ Database Fixtures ##################
@pytest.fixture()
def database_test_directory(test_data_directory):
    return os.path.join(test_data_directory, "DatabaseTests")


@pytest.fixture()
def database_chemsmart_file(database_test_directory):
    return os.path.join(database_test_directory, "chemsmart.db")


@pytest.fixture()
def database_chemsmart_xtb_file(database_test_directory):
    return os.path.join(database_test_directory, "chemsmart_xtb.db")


@pytest.fixture()
def database_ase_file(database_test_directory):
    return os.path.join(database_test_directory, "ase.db")


@pytest.fixture()
def database_empty_file(database_test_directory):
    return os.path.join(database_test_directory, "empty.db")
