from pathlib import Path

from click.testing import CliRunner

from chemsmart.cli.run import run


def _write_signature_file(path: Path, program: str):
    signatures = {
        "gaussian": "Gaussian, Inc.\n",
        "orca": "* O   R   C   A *\n",
        "unknown": "Some random text\n",
    }
    path.write_text(signatures[program])


def _build_outputs(tmp_path: Path, program: str):
    names = [
        "ha.log",
        "a.log",
        "hb.log",
        "b.log",
        "has.log",
        "as.log",
        "hbs.log",
        "bs.log",
    ]
    files = {}
    for n in names:
        p = tmp_path / n
        _write_signature_file(p, program)
        files[n] = str(p)
    return files


def _invoke_pka(runner, files):
    return runner.invoke(
        run,
        [
            "pka",
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-hr",
            files["hb.log"],
            "-r",
            files["b.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
            "--href-solv",
            files["hbs.log"],
            "--ref-solv",
            files["bs.log"],
            "-rp",
            "6.75",
        ],
    )


def test_run_pka_detects_gaussian_and_dispatches(tmp_path, monkeypatch):
    files = _build_outputs(tmp_path, "gaussian")
    called = {}

    def _fake_print(**kwargs):
        called["kwargs"] = kwargs

    from chemsmart.io.gaussian.output import Gaussian16pKaOutput

    monkeypatch.setattr(Gaussian16pKaOutput, "print_pka_summary", _fake_print)

    runner = CliRunner()
    result = _invoke_pka(runner, files)

    assert result.exit_code == 0
    assert "kwargs" in called
    assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
    assert called["kwargs"]["a_solv_file"] == files["as.log"]
    assert called["kwargs"]["pka_reference"] == 6.75


def test_run_pka_detects_orca_and_dispatches(tmp_path, monkeypatch):
    files = _build_outputs(tmp_path, "orca")
    called = {}

    def _fake_print(**kwargs):
        called["kwargs"] = kwargs

    from chemsmart.io.orca.output import ORCApKaOutput

    monkeypatch.setattr(ORCApKaOutput, "print_pka_summary", _fake_print)

    runner = CliRunner()
    result = _invoke_pka(runner, files)

    assert result.exit_code == 0
    assert "kwargs" in called
    assert called["kwargs"]["hb_gas_file"] == files["hb.log"]


def test_run_pka_mixed_programs_raises(tmp_path):
    files = _build_outputs(tmp_path, "gaussian")
    _write_signature_file(Path(files["bs.log"]), "orca")

    runner = CliRunner()
    result = _invoke_pka(runner, files)

    assert result.exit_code != 0
    assert "mixed program types" in result.output


def test_run_pka_unknown_program_raises(tmp_path):
    files = _build_outputs(tmp_path, "unknown")

    runner = CliRunner()
    result = _invoke_pka(runner, files)

    assert result.exit_code != 0
    assert "Could not detect output-file program type" in result.output
