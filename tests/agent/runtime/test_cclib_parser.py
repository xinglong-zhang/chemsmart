from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.runtime import cclib_parser
from chemsmart.agent.runtime.calculations import inspect_output

REPO_ROOT = Path(__file__).parents[3]


def test_cclib_parses_gaussian_ts_and_preserves_native_termination():
    output = (
        REPO_ROOT
        / "tests"
        / "data"
        / "GaussianTests"
        / "outputs"
        / "Pd_insertion_ts_r.log"
    )

    summary = inspect_output(output, program="gaussian", kind="ts")

    assert summary["parser_backend"] == "native+cclib"
    assert summary["normal_termination"] is True
    assert summary["atom_count"] == 42
    assert summary["imaginary_frequency_count"] == 1
    assert summary["parsed_charge"] == 0
    assert summary["parsed_multiplicity"] == 1


def test_cclib_parses_orca_optimization_fixture():
    output = (
        REPO_ROOT
        / "tests"
        / "data"
        / "ORCATests"
        / "outputs"
        / "water_opt.out"
    )

    summary = inspect_output(output, program="orca", kind="opt")

    assert summary["parser_backend"] == "native+cclib"
    assert summary["normal_termination"] is True
    assert summary["atom_count"] == 3
    assert summary["frequency_count"] == 3


def test_cclib_parser_failure_is_evidence_not_calculation_failure(
    tmp_path, monkeypatch
):
    output = tmp_path / "orca.out"
    output.write_text(
        "FINAL SINGLE POINT ENERGY -76.269459371830\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )

    class BrokenParser:
        def parse(self):
            raise IndexError("ORCA 6 SCF layout")

    monkeypatch.setattr(
        cclib_parser, "ccopen", lambda *args, **kwargs: BrokenParser()
    )

    summary = inspect_output(output, program="orca", kind="sp")

    assert summary["normal_termination"] is True
    assert summary["energy"] == pytest.approx(-76.269459371830)
    assert summary["parser_backend"] == "native"
    assert summary["parser_warnings"] == [
        "cclib IndexError: ORCA 6 SCF layout"
    ]


def test_missing_cclib_keeps_native_result_parser(tmp_path, monkeypatch):
    output = tmp_path / "orca.out"
    output.write_text(
        "FINAL SINGLE POINT ENERGY -76.269459371830\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(cclib_parser, "ccopen", None)
    monkeypatch.setattr(cclib_parser, "convertor", None)

    summary = inspect_output(output, program="orca", kind="sp")

    assert summary["normal_termination"] is True
    assert summary["energy"] == pytest.approx(-76.269459371830)
    assert summary["parser_backend"] == "native"
    assert summary["parser_warnings"] == [
        "cclib is not installed; using the built-in result parser. Install "
        "chemsmart[result-parsers] for normalized advanced properties."
    ]
