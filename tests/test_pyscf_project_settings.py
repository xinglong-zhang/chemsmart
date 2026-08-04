"""Compute-free contracts for canonical PySCF project stages."""

from pathlib import Path

import pytest

from chemsmart.jobs.pyscf.opt import PySCFOptJob
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.settings.pyscf import PySCFProjectSettings
from chemsmart.settings.pyscf import materialize_canonical_project_sections


def _project(tmp_path, body):
    path = Path(tmp_path) / "project.yaml"
    path.write_text(body, encoding="utf-8")
    return PySCFProjectSettings.from_yaml(path)


def test_canonical_opt_cannot_hide_an_in_process_hessian(tmp_path):
    project = _project(
        tmp_path,
        "opt:\n  functional: b3lyp\n  basis: def2-svp\n  freq: true\n",
    )

    with pytest.raises(ValueError, match="separate hess stage"):
        project.opt_settings()


def test_legacy_opt_freq_migrates_to_explicit_opt_and_hess_nodes(tmp_path):
    project = _project(
        tmp_path,
        "gas:\n  functional: b3lyp\n  basis: def2-svp\n  freq: true\n",
    )

    assert project.opt_settings().freq is False
    assert project.hess_settings().freq is True
    assert PySCFOptJob.stages.fget(None) == ["scf", "opt"]
    canonical = project.render_canonical_yaml(jobtypes=("opt", "hess"))
    assert "gas:" not in canonical
    assert "opt:" in canonical
    assert "hess:" in canonical


def test_canonical_sections_make_loader_effective_defaults_explicit(tmp_path):
    project = _project(
        tmp_path,
        "sp:\n  functional: b3lyp\n  basis: def2-svp\n"
        "opt:\n  functional: b3lyp\n  basis: def2-svp\n"
        "hess:\n  functional: b3lyp\n  basis: def2-svp\n",
    )

    sections = project.canonical_sections()

    for stage in ("sp", "opt", "hess"):
        assert sections[stage]["functional"] == "b3lyp"
        assert sections[stage]["basis"] == "def2-svp"
        assert sections[stage]["density_fit"] is False
        assert sections[stage]["engine"] == "cpu"
        assert sections[stage]["scf_tol"] is None
        assert sections[stage]["scf_maxiter"] is None
    assert sections["sp"]["freq"] is False
    assert sections["opt"]["freq"] is False
    assert sections["hess"]["freq"] is True
    assert sections["opt"]["opt_solver"] == "geometric"
    assert sections["opt"]["opt_maxsteps"] == 100
    assert "opt_solver" not in sections["sp"]


def test_canonical_materialization_preserves_explicit_values_and_legacy_input():
    sections = materialize_canonical_project_sections(
        {
            "solv": {
                "functional": "pbe0",
                "basis": "def2-tzvp",
                "density_fit": True,
                "scf_tol": 1.0e-10,
                "scf_maxiter": 75,
            },
            "gas": {
                "functional": "pbe0",
                "basis": "def2-tzvp",
                "opt_solver": "berny",
                "opt_maxsteps": 60,
            },
        }
    )

    assert tuple(sections) == ("sp", "opt", "hess")
    assert sections["sp"]["density_fit"] is True
    assert sections["sp"]["scf_tol"] == 1.0e-10
    assert sections["sp"]["scf_maxiter"] == 75
    assert sections["opt"]["opt_solver"] == "berny"
    assert sections["opt"]["opt_maxsteps"] == 60
    assert sections["hess"]["freq"] is True


def test_td_project_materializes_explicit_bounded_response_fields(tmp_path):
    project = _project(
        tmp_path,
        "td:\n"
        "  functional: b3lyp\n"
        "  basis: def2-svp\n"
        "  response_method: tddft\n"
        "  state_manifold: singlet\n"
        "  nstates: 6\n",
    )

    settings = project.td_settings()

    assert settings.jobtype == "td"
    assert settings.response_method == "tddft"
    assert settings.state_manifold == "singlet"
    assert settings.nstates == 6
    assert settings.freq is False


def test_td_settings_reject_gpu_response_preview():
    settings = PySCFJobSettings(
        jobtype="td",
        functional="b3lyp",
        basis="def2-svp",
        response_method="tda",
        state_manifold="singlet",
        nstates=3,
        charge=0,
        multiplicity=1,
        engine="gpu",
    )

    with pytest.raises(ValueError, match="CPU-only preview"):
        settings.validate()
