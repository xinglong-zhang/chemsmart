import importlib
from pathlib import Path
from typing import cast

from click.core import BaseCommand
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub


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

    def _fake_print(*args, **kwargs):
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

    def _fake_print(*args, **kwargs):
        called["kwargs"] = kwargs

    from chemsmart.io.orca.output import ORCApKaOutput

    monkeypatch.setattr(ORCApKaOutput, "print_pka_summary", _fake_print)

    runner = CliRunner()
    result = _invoke_pka(runner, files)

    assert result.exit_code == 0
    assert "kwargs" in called
    assert called["kwargs"]["href_gas_file"] == files["hb.log"]


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


def test_sub_orca_pka_batch_array_reconstructs_per_job_cli_args(
    tmp_path, monkeypatch
):
    orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.settings.server import Server

    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    table = tmp_path / "batch.xyz"
    table.write_text(
        "filepath proton_index charge multiplicity\n"
        f"{acid1} 2 0 1\n"
        f"{acid2} 2 1 2\n"
    )

    config_root = tmp_path / "chemsmart_cfg"
    orca_cfg_dir = config_root / "orca"
    orca_cfg_dir.mkdir(parents=True)
    (orca_cfg_dir / "test.yaml").write_text(
        "gas:\n"
        "  functional: B3LYP\n"
        "  basis: def2-SVP\n"
        "solv:\n"
        "  functional: B3LYP\n"
        "  basis: def2-SVP\n"
        "  freq: false\n"
        "  solvent_model: smd\n"
        "  solvent_id: water\n"
    )
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

    captured = {}

    fake_server = Server(name="dummy")
    real_from_filepath = Molecule.from_filepath

    def _fake_from_filepath(filepath, *args, **kwargs):
        if str(filepath) == str(table):
            placeholder = Molecule(
                symbols=["C", "H"],
                positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                charge=0,
                multiplicity=1,
            )
            if kwargs.get("return_list"):
                return [placeholder]
            return placeholder
        return real_from_filepath(filepath, *args, **kwargs)

    def _fake_submit_array_job(
        jobs, num_nodes=None, test=False, cli_args=None, **kwargs
    ):
        captured["jobs"] = jobs
        captured["num_nodes"] = num_nodes
        captured["test"] = test
        captured["cli_args"] = cli_args

    monkeypatch.setattr(
        fake_server, "submit_array_job", _fake_submit_array_job
    )
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    monkeypatch.setattr(
        orca_cli.Molecule, "from_filepath", _fake_from_filepath
    )

    runner = CliRunner()
    result = runner.invoke(
        cast(BaseCommand, sub),
        [
            "--server",
            "dummy",
            "--num-nodes",
            "2",
            "--test",
            "orca",
            "--project",
            "test",
            "--filename",
            str(table),
            "pka",
            "--scheme",
            "direct",
            "batch",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["num_nodes"] == 2
    assert captured["test"] is True
    assert len(captured["jobs"]) == 2
    assert len(captured["cli_args"]) == 2

    first_args, second_args = captured["cli_args"]

    assert isinstance(first_args, list)
    assert isinstance(second_args, list)
    assert first_args != second_args

    assert str(acid1) in first_args
    assert str(acid2) in second_args
    assert str(acid2) not in first_args
    assert str(acid1) not in second_args

    assert "submit" in first_args
    assert "submit" in second_args
    assert "batch" not in first_args
    assert "batch" not in second_args

    assert "--charge" in first_args
    assert first_args[first_args.index("--charge") + 1] == "0"
    assert "--multiplicity" in first_args
    assert first_args[first_args.index("--multiplicity") + 1] == "1"
    assert "--proton-index" in first_args
    assert first_args[first_args.index("--proton-index") + 1] == "2"

    assert "--charge" in second_args
    assert second_args[second_args.index("--charge") + 1] == "1"
    assert "--multiplicity" in second_args
    assert second_args[second_args.index("--multiplicity") + 1] == "2"
    assert "--proton-index" in second_args
    assert second_args[second_args.index("--proton-index") + 1] == "2"
