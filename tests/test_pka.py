import importlib
from pathlib import Path

import click
import pytest
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


def _invoke_pka_direct(runner, files, delta_g_proton=-265.9):
    return runner.invoke(
        run,
        [
            "pka",
            "-s",
            "direct",
            "-dG",
            str(delta_g_proton),
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
        ],
    )


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


def _require_backend_pka_subcommand(command_group, backend):
    runner = CliRunner()
    result = runner.invoke(command_group, [backend, "--help"])
    assert result.exit_code == 0, result.output
    if "\n  pka" not in result.output:
        pytest.skip(
            f"{backend} backend pka subcommand is not registered in this build."
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


def test_run_pka_direct_analyze_dispatches(tmp_path, monkeypatch):
    files = _build_outputs(tmp_path, "gaussian")
    called = {}

    def _fake_print(*args, **kwargs):
        called["kwargs"] = kwargs

    from chemsmart.io.gaussian.output import Gaussian16pKaOutput

    monkeypatch.setattr(Gaussian16pKaOutput, "print_pka_summary", _fake_print)

    runner = CliRunner()
    result = _invoke_pka_direct(runner, files, delta_g_proton=-270.0)

    assert result.exit_code == 0
    assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
    assert called["kwargs"]["delta_G_proton"] == -270.0
    assert called["kwargs"]["scheme"] == "direct"


def test_run_pka_direct_requires_delta_g_proton(tmp_path):
    files = _build_outputs(tmp_path, "gaussian")
    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            "pka",
            "-s",
            "direct",
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
        ],
    )
    assert result.exit_code != 0
    assert "-dG/--delta-g-proton is required" in result.output


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
    assert isinstance(result.exception, ValueError)
    assert "mixed program types" in str(result.exception).lower()


def test_run_pka_unknown_program_raises(tmp_path):
    files = _build_outputs(tmp_path, "unknown")

    runner = CliRunner()
    result = _invoke_pka(runner, files)

    assert result.exit_code != 0
    assert isinstance(result.exception, ValueError)
    assert (
        "could not detect output-file program type"
        in str(result.exception).lower()
    )


def test_sub_orca_pka_batch_reconstructs_per_job_cli_args(
    tmp_path, monkeypatch
):
    _require_backend_pka_subcommand(sub, "orca")
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

    captured = {"submissions": []}

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

    def _fake_submit(job, test=False, cli_args=None, **kwargs):
        captured["submissions"].append((job, test, cli_args))

    monkeypatch.setattr(fake_server, "submit", _fake_submit)
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    monkeypatch.setattr(
        orca_cli.Molecule, "from_filepath", _fake_from_filepath
    )

    runner = CliRunner()
    result = runner.invoke(
        sub,
        [
            "--server",
            "dummy",
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
    assert len(captured["submissions"]) == 2
    first_job, first_test, first_args = captured["submissions"][0]
    second_job, second_test, second_args = captured["submissions"][1]
    assert first_test is True
    assert second_test is True
    assert isinstance(first_args, list)
    assert isinstance(second_args, list)
    # Per-entry submit scripts should execute a single-row submission.
    assert "submit" in first_args
    assert "batch" not in first_args
    assert str(table) not in first_args


def test_sub_orca_pka_batch_rewrites_per_entry_file_args(
    tmp_path, monkeypatch
):
    _require_backend_pka_subcommand(sub, "orca")
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

    captured = {"submissions": []}

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

    def _fake_submit(job, test=False, cli_args=None, **kwargs):
        captured["submissions"].append((job, test, cli_args))

    monkeypatch.setattr(fake_server, "submit", _fake_submit)
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    monkeypatch.setattr(
        orca_cli.Molecule, "from_filepath", _fake_from_filepath
    )

    runner = CliRunner()
    result = runner.invoke(
        sub,
        [
            "--server",
            "dummy",
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
    assert len(captured["submissions"]) == 2

    first_args = captured["submissions"][0][2]
    second_args = captured["submissions"][1][2]

    assert str(table) not in first_args
    assert str(table) not in second_args
    assert str(acid1) in first_args
    assert str(acid2) in second_args
    # Rewritten row-level options must be placed in the correct command scope:
    # --charge/--multiplicity on backend command and --proton-index on pka.
    assert "--charge" in first_args
    assert "--multiplicity" in first_args
    assert "--proton-index" in first_args
    assert "submit" in first_args
    assert "batch" not in first_args
    assert first_args.index("--charge") < first_args.index("pka")
    assert first_args.index("--multiplicity") < first_args.index("pka")
    assert first_args.index("--proton-index") < first_args.index("submit")


def test_sub_orca_pka_batch_shared_reference_loaded_once(
    tmp_path, monkeypatch
):
    _require_backend_pka_subcommand(sub, "orca")
    orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.orca.settings import ORCApKaJobSettings
    from chemsmart.settings.server import Server

    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    reference = tmp_path / "ref.xyz"
    reference.write_text("2\nref\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

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

    captured = {"submissions": []}
    reference_pair_call_count = {"count": 0}

    fake_server = Server(name="dummy")
    real_from_filepath = Molecule.from_filepath
    real_reference_pair_molecules = ORCApKaJobSettings.reference_pair_molecules

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

    def _counting_reference_pair(self):
        reference_pair_call_count["count"] += 1
        return real_reference_pair_molecules(self)

    def _fake_submit(job, test=False, cli_args=None, **kwargs):
        captured["submissions"].append((job, test, cli_args))

    monkeypatch.setattr(fake_server, "submit", _fake_submit)
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    monkeypatch.setattr(
        orca_cli.Molecule, "from_filepath", _fake_from_filepath
    )
    monkeypatch.setattr(
        ORCApKaJobSettings,
        "reference_pair_molecules",
        _counting_reference_pair,
    )

    runner = CliRunner()
    result = runner.invoke(
        sub,
        [
            "--server",
            "dummy",
            "--test",
            "orca",
            "--project",
            "test",
            "--filename",
            str(table),
            "pka",
            "--scheme",
            "proton exchange",
            "--reference",
            str(reference),
            "--reference-proton-index",
            "2",
            "--reference-charge",
            "0",
            "--reference-multiplicity",
            "1",
            "batch",
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(captured["submissions"]) == 2
    # In "sub ... pka batch", jobs are not executed; only submission scripts are
    # generated, so reference molecules are not built at this stage.
    assert reference_pair_call_count["count"] == 0


def test_sub_orca_pka_batch_first_exchange_rest_direct(tmp_path, monkeypatch):
    _require_backend_pka_subcommand(sub, "orca")
    orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.settings.server import Server

    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    reference = tmp_path / "ref.xyz"
    reference.write_text("2\nref\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

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

    captured = {"submissions": []}

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

    def _fake_submit(job, test=False, cli_args=None, **kwargs):
        captured["submissions"].append((job, test, cli_args))

    monkeypatch.setattr(fake_server, "submit", _fake_submit)
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    monkeypatch.setattr(
        orca_cli.Molecule, "from_filepath", _fake_from_filepath
    )

    runner = CliRunner()
    result = runner.invoke(
        sub,
        [
            "--server",
            "dummy",
            "--test",
            "orca",
            "--project",
            "test",
            "--filename",
            str(table),
            "pka",
            "--scheme",
            "proton exchange",
            "--reference",
            str(reference),
            "--reference-proton-index",
            "2",
            "--reference-charge",
            "0",
            "--reference-multiplicity",
            "1",
            "batch",
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(captured["submissions"]) == 2

    first_job, _, first_args = captured["submissions"][0]
    second_job, _, second_args = captured["submissions"][1]

    assert first_job.settings.scheme == "proton exchange"
    assert second_job.settings.scheme == "direct"

    assert "--scheme" in first_args
    assert first_args[first_args.index("--scheme") + 1] == "proton exchange"
    assert "--reference" in first_args

    assert "--scheme" in second_args
    assert second_args[second_args.index("--scheme") + 1] == "direct"
    assert "--reference" not in second_args
    assert "--reference-proton-index" not in second_args
    assert "--reference-charge" not in second_args
    assert "--reference-multiplicity" not in second_args


def test_run_gaussian_pka_help_is_submission_only(
    single_molecule_xyz_file,
):
    _require_backend_pka_subcommand(run, "gaussian")
    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            "--no-scratch",
            "--fake",
            "gaussian",
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "pka",
            "--help",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "\n  submit" in result.output
    assert "\n  batch" in result.output
    assert "\n  analyze" not in result.output
    assert "\n  thermo" not in result.output
    assert "\n  batch-analyze" not in result.output


def test_run_orca_pka_help_is_submission_only(
    single_molecule_xyz_file,
):
    _require_backend_pka_subcommand(run, "orca")
    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            "--no-scratch",
            "--fake",
            "orca",
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "pka",
            "--help",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "\n  submit" in result.output
    assert "\n  batch" in result.output
    assert "\n  analyze" not in result.output
    assert "\n  thermo" not in result.output
    assert "\n  batch-analyze" not in result.output


def test_run_pka_help_keeps_output_analysis_commands():
    runner = CliRunner()
    result = runner.invoke(run, ["--no-scratch", "--fake", "pka", "--help"])

    assert result.exit_code == 0, result.output
    assert "analyze" in result.output
    assert "batch-analyze" in result.output


def _write_test_backend_project(tmp_path, backend):
    config_root = tmp_path / "chemsmart_cfg"
    backend_cfg_dir = config_root / backend
    backend_cfg_dir.mkdir(parents=True)
    (backend_cfg_dir / "test.yaml").write_text(
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
    return config_root


def _build_pka_batch_table(tmp_path):
    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    table = tmp_path / "pka_scale.csv"
    table.write_text(
        "filepath,proton_index,charge,multiplicity\n"
        f"{acid1},2,0,1\n"
        f"{acid2},2,1,2\n"
    )
    return table


@pytest.mark.parametrize("backend", ["gaussian", "orca"])
def test_run_pka_batch_table_processing(tmp_path, monkeypatch, backend):
    """pKa table batch returns multiple jobs; run executes each locally."""
    _require_backend_pka_subcommand(run, backend)
    table = _build_pka_batch_table(tmp_path)
    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

    captured = {"runs": []}

    from chemsmart.jobs.job import Job

    def _fake_run(self):
        captured["runs"].append(self.label)

    monkeypatch.setattr(Job, "run", _fake_run)

    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            "--no-scratch",
            "--fake",
            backend,
            "-p",
            "test",
            "-f",
            str(table),
            "pka",
            "-s",
            "direct",
            "batch",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Batch job submission is not supported" not in result.output
    assert len(captured["runs"]) == 2
    if backend == "gaussian":
        assert set(captured["runs"]) == {"acid1", "acid2"}
    else:
        assert set(captured["runs"]) == {"acid1_pka", "acid2_pka"}


@pytest.mark.parametrize("backend", ["gaussian", "orca"])
def test_run_pka_batch_defaults_to_no_scratch(tmp_path, monkeypatch, backend):
    """run should not require a scratch directory when --scratch is omitted."""
    _require_backend_pka_subcommand(run, backend)
    table = _build_pka_batch_table(tmp_path)
    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

    missing_scratch = tmp_path / "missing_scratch"
    from chemsmart.jobs import runner as runner_module

    monkeypatch.setattr(
        runner_module.user_settings, "scratch", str(missing_scratch)
    )

    captured = {"runs": []}

    from chemsmart.jobs.job import Job

    def _fake_run(self):
        captured["runs"].append(self.label)

    monkeypatch.setattr(Job, "run", _fake_run)

    runner = CliRunner()
    result = runner.invoke(
        run,
        [
            "--fake",
            backend,
            "-p",
            "test",
            "-f",
            str(table),
            "pka",
            "-s",
            "direct",
            "batch",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Specified scratch dir does not exist" not in result.output
    assert len(captured["runs"]) == 2


def test_run_rejects_non_job_batch_payload(pbs_server):
    """Scheduler-style batch payloads that are not Job lists stay blocked."""

    from chemsmart.cli.run import process_pipeline
    from chemsmart.jobs.runner import JobRunner

    ctx = click.Context(run)
    ctx.ensure_object(dict)
    ctx.obj["jobrunner"] = JobRunner(server=pbs_server, fake=True)

    with pytest.raises(
        ValueError, match="Batch job submission is not supported"
    ):
        process_pipeline.__wrapped__(ctx, ["not-a-job", "also-not-a-job"])
