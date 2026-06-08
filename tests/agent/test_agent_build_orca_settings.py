from __future__ import annotations

from chemsmart.agent.tools import build_orca_settings


def test_dlpno_ccsd_t_routes_to_ab_initio_not_functional():
    settings = build_orca_settings("DLPNO-CCSD(T)", "def2-TZVP")

    assert settings.ab_initio == "DLPNO-CCSD(T)"
    assert settings.functional is None
    assert settings.basis == "def2-TZVP"


def test_mp2_routes_to_ab_initio():
    settings = build_orca_settings("MP2", "def2-SVP")

    assert settings.ab_initio == "MP2"
    assert settings.functional is None
    assert settings.basis == "def2-SVP"


def test_dft_functional_not_touched():
    settings = build_orca_settings("B3LYP", "def2-SVP")

    assert settings.functional == "B3LYP"
    assert settings.ab_initio is None
    assert settings.basis == "def2-SVP"
