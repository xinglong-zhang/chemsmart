"""The shared ``monkeypatch`` fixture must isolate ``~`` on every platform.

Tests redirect the ChemSmart user tree with ``setenv("HOME", tmp_path)``.
Windows resolves ``~`` through ``USERPROFILE``, so without mirroring, those
tests wrote into the real profile: the wizard suite failed with
``FileExistsError`` on leftovers from earlier runs instead of exercising a
clean sandbox.
"""

from __future__ import annotations

import os
from pathlib import Path


def test_setenv_home_also_redirects_expanduser(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))

    assert Path.home() == tmp_path
    assert Path(os.path.expanduser("~")) == tmp_path
    assert os.environ["USERPROFILE"] == str(tmp_path)


def test_home_isolation_is_undone_after_the_test(monkeypatch, tmp_path):
    original_home = os.environ.get("HOME")
    original_profile = os.environ.get("USERPROFILE")

    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.undo()

    assert os.environ.get("HOME") == original_home
    assert os.environ.get("USERPROFILE") == original_profile


def test_other_variables_are_untouched(monkeypatch, tmp_path):
    monkeypatch.setenv("CHEMSMART_TEST_SENTINEL", "value")

    assert os.environ["CHEMSMART_TEST_SENTINEL"] == "value"
    assert os.environ.get("USERPROFILE") != "value"


def test_delenv_home_clears_the_alias(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("HOME")

    assert "HOME" not in os.environ
    assert "USERPROFILE" not in os.environ
