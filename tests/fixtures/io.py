"""IO fixtures."""

import pytest


############ IO Fixtures ####################################
@pytest.fixture
def temporary_working_dir(tmp_path, monkeypatch):
    """Run a test in an isolated temporary working directory."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture()
def captured():
    """Return an empty dictionary for capturing mock call arguments."""
    return {}
