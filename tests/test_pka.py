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


class _FakeThermochemistry:
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.electronic_energy = -627.0
        self.qrrho_gibbs_free_energy = -628.0


def _install_fake_thermochemistry(monkeypatch, constructed=None):
    constructed = [] if constructed is None else constructed

    class _TrackingFakeThermochemistry(_FakeThermochemistry):
        def __init__(self, filename, **kwargs):
            constructed.append(Path(filename).name)
            super().__init__(filename, **kwargs)

    monkeypatch.setattr(
        "chemsmart.cli.pka.Thermochemistry",
        _TrackingFakeThermochemistry,
    )
    return constructed


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


def _setup_sub_pka_batch_test(tmp_path, monkeypatch, backend):
    """Shared fixtures for sub ... pka batch submission tests."""
    monkeypatch.chdir(tmp_path)
    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    table = tmp_path / "pka_scale.csv"
    table.write_text(
        "structure,filepath,proton_index,charge,multiplicity\n"
        f"acid1,{acid1},2,0,1\n"
        f"acid2,{acid2},2,1,2\n"
    )

    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

    from chemsmart.settings.server import Server

    fake_server = Server(name="dummy")
    captured = {"submissions": []}

    def _fake_submit_array(
        jobs,
        array_concurrency=None,
        test=False,
        cli_args=None,
        batch_label=None,
        **kw,
    ):
        captured["submissions"].append(
            {
                "jobs": list(jobs),
                "test": test,
                "cli_args": cli_args,
                "batch_label": batch_label,
                "array_concurrency": array_concurrency,
            }
        )

    fake_server.submit_array_job = _fake_submit_array
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    return table, captured


def _submitted_jobs(captured):
    """Return job objects recorded by fake submit handlers."""
    jobs = []
    for entry in captured["submissions"]:
        if isinstance(entry, dict):
            jobs.extend(entry["jobs"])
        else:
            jobs.append(entry[0])
    return jobs


def _submitted_pka_child_jobs(captured):
    """Unwrap a submitted pKa batch container into its child jobs."""
    from chemsmart.jobs.batch import BatchJob

    if not captured["submissions"]:
        return []
    first = captured["submissions"][0]
    if isinstance(first, dict):
        return list(first["jobs"])
    job = first[0]
    if isinstance(job, BatchJob):
        return list(job.jobs)
    if isinstance(job, list):
        return job
    return [job]


def _submitted_cli_args(captured):
    """Return cli_args from the first captured array/single submission."""
    if not captured["submissions"]:
        return None
    first = captured["submissions"][0]
    if isinstance(first, dict):
        return first["cli_args"]
    return first[2]


def _attach_fake_array_server(monkeypatch, tmp_path, captured=None):
    """Patch Server.from_servername to capture submit_array_job calls."""
    from chemsmart.settings.server import Server

    monkeypatch.chdir(tmp_path)
    if captured is None:
        captured = {"submissions": []}
    fake_server = Server(name="dummy")

    def _fake_submit_array(
        jobs,
        array_concurrency=None,
        test=False,
        cli_args=None,
        batch_label=None,
        **kw,
    ):
        captured["submissions"].append(
            {
                "jobs": list(jobs),
                "test": test,
                "cli_args": cli_args,
                "batch_label": batch_label,
                "array_concurrency": array_concurrency,
            }
        )

    def _fake_submit(job, test=False, cli_args=None, **kw):
        captured["submissions"].append(
            {
                "jobs": [job],
                "test": test,
                "cli_args": cli_args,
                "batch_label": None,
                "array_concurrency": None,
            }
        )

    fake_server.submit_array_job = _fake_submit_array
    fake_server.submit = _fake_submit
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    return fake_server, captured


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


def test_rewrite_pka_batch_cli_args_replaces_table_with_submit_row():
    from chemsmart.cli.pka import rewrite_pka_batch_cli_args

    shared = [
        "gaussian",
        "-f",
        "table.csv",
        "pka",
        "--scheme",
        "proton exchange",
        "--reference",
        "ref.xyz",
        "batch",
    ]
    entry = {
        "filepath": "acid2.xyz",
        "proton_index": 3,
        "charge": 1,
        "multiplicity": 2,
        "scheme": "direct",
        "label": "acid2",
    }
    rewritten = rewrite_pka_batch_cli_args(shared, entry)
    assert "batch" not in rewritten
    assert "submit" in rewritten
    assert "acid2.xyz" in rewritten
    assert "table.csv" not in rewritten
    assert rewritten[rewritten.index("--proton-index") + 1] == "3"
    assert rewritten[rewritten.index("--charge") + 1] == "1"
    assert rewritten[rewritten.index("--scheme") + 1] == "direct"
    assert "--reference" not in rewritten


class TestPKa:
    """pKa CLI, batch submission, and job workflow tests."""

    def test_run_pka_detects_gaussian_and_dispatches(
        self, tmp_path, monkeypatch
    ):
        files = _build_outputs(tmp_path, "gaussian")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code == 0
        assert "kwargs" in called
        assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
        assert called["kwargs"]["a_solv_file"] == files["as.log"]
        assert called["kwargs"]["pka_reference"] == 6.75

    def test_run_pka_direct_analyze_dispatches(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka_direct(runner, files, delta_g_proton=-270.0)

        assert result.exit_code == 0
        assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
        assert called["kwargs"]["delta_G_proton"] == -270.0
        assert called["kwargs"]["scheme"] == "direct"

    def test_run_pka_direct_requires_delta_g_proton(self, tmp_path):
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

    def test_run_pka_detects_orca_and_dispatches(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "orca")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code == 0
        assert "kwargs" in called
        assert called["kwargs"]["href_gas_file"] == files["hb.log"]

    def test_run_pka_mixed_programs_analyze(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        _write_signature_file(Path(files["bs.log"]), "orca")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code == 0, result.output
        assert called["kwargs"]["href_solv_file"] == files["hbs.log"]
        assert called["kwargs"]["ref_solv_file"] == files["bs.log"]

    def test_run_pka_batch_analyze_orca_outputs(self, tmp_path, monkeypatch):
        """batch-analyze should build Thermochemistry objects for ORCA files."""
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.out").write_text(
                "* O   R   C   A *\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.out").write_text("* O   R   C   A *\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.out,ref_A_opt.out,ref_HA_sp.out,ref_A_sp.out,10.6\n"
        )

        constructed = _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "-T",
                "333.15",
                "-csg",
                "100",
                "-ch",
                "100",
                "batch-analyze",
                "-o",
                str(table),
            ],
        )

        assert result.exit_code == 0, result.output
        assert "target_pka_HA_sp.out" in constructed
        assert "ref_HA_sp.out" in constructed
        assert "pKa" in result.output

    def test_run_pka_batch_analyze_mixed_gaussian_orca(
        self, tmp_path, monkeypatch
    ):
        """batch-analyze stays program-agnostic via Thermochemistry(filename=...)."""
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.out").write_text(
                "* O   R   C   A *\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.log").write_text("Gaussian, Inc.\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.log,ref_A_opt.log,ref_HA_sp.log,ref_A_sp.log,10.6\n"
        )

        constructed = _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "batch-analyze",
                "-o",
                str(table),
            ],
        )

        assert result.exit_code == 0, result.output
        assert "target_pka_HA_sp.out" in constructed
        assert "ref_HA_sp.log" in constructed
        assert "pKa" in result.output

    def test_run_pka_batch_analyze_gaussian_outputs(
        self, tmp_path, monkeypatch
    ):
        """batch-analyze should preserve Gaussian table behavior."""
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.log").write_text(
                "Gaussian, Inc.\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.log").write_text("Gaussian, Inc.\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.log,ref_A_opt.log,ref_HA_sp.log,ref_A_sp.log,6.75\n"
        )

        constructed = _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            ["pka", "batch-analyze", "-o", str(table)],
        )

        assert result.exit_code == 0, result.output
        assert "target_pka_HA_opt.log" in constructed
        assert "ref_HA_sp.log" in constructed
        assert "pKa" in result.output

    def test_pka_thermochemistry_missing_scf_energy(
        self, tmp_path, monkeypatch
    ):
        class _MissingScfThermochemistry:
            electronic_energy = None
            qrrho_gibbs_free_energy = -1.0

            def __init__(self, filename, **kwargs):
                pass

        monkeypatch.setattr(
            "chemsmart.cli.pka.Thermochemistry",
            _MissingScfThermochemistry,
        )

        from chemsmart.cli.pka import pka_solvent_scf_energy

        with pytest.raises(ValueError, match="Could not extract SCF energy"):
            pka_solvent_scf_energy(str(tmp_path / "missing.out"))

    def test_pka_thermochemistry_missing_qh_gibbs(self, tmp_path, monkeypatch):
        class _MissingQhThermochemistry:
            electronic_energy = -1.0
            qrrho_gibbs_free_energy = None

            def __init__(self, filename, **kwargs):
                pass

        monkeypatch.setattr(
            "chemsmart.cli.pka.Thermochemistry",
            _MissingQhThermochemistry,
        )

        from chemsmart.cli.pka import pka_gas_phase_data

        with pytest.raises(
            ValueError,
            match="Could not extract quasi-harmonic Gibbs free energy",
        ):
            pka_gas_phase_data(str(tmp_path / "gas.out"))

    def test_run_pka_unparseable_output_raises(self, tmp_path):
        """analyze no longer pre-detects program type; parsing fails on bad files."""
        files = _build_outputs(tmp_path, "unknown")

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code != 0

    def test_sub_orca_pka_batch_reconstructs_per_job_cli_args(
        self, tmp_path, monkeypatch
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule

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

        fake_server, captured = _attach_fake_array_server(
            monkeypatch, tmp_path, captured
        )
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
        assert len(captured["submissions"]) == 1
        submission = captured["submissions"][0]
        assert submission["test"] is True
        cli_args = submission["cli_args"]
        assert isinstance(cli_args, list)
        assert len(cli_args) == 2
        assert all("submit" in args for args in cli_args)
        assert str(acid1) in cli_args[0]
        assert str(acid2) in cli_args[1]
        children = _submitted_pka_child_jobs(captured)
        assert len(children) == 2

    def test_sub_orca_pka_batch_rewrites_per_entry_file_args(
        self, tmp_path, monkeypatch
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule

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
        _attach_fake_array_server(monkeypatch, tmp_path, captured)
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
        assert len(captured["submissions"]) == 1
        children = _submitted_pka_child_jobs(captured)
        assert len(children) == 2
        cli_args = _submitted_cli_args(captured)
        assert len(cli_args) == 2
        assert str(acid1) in cli_args[0]
        assert str(acid2) in cli_args[1]
        assert "submit" in cli_args[0]
        assert "batch" not in cli_args[0]
        assert str(table) not in cli_args[0]

    def test_sub_orca_pka_batch_shared_reference_loaded_once(
        self, tmp_path, monkeypatch
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

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

        _attach_fake_array_server(monkeypatch, tmp_path, captured)
        real_from_filepath = Molecule.from_filepath
        real_reference_pair_molecules = (
            ORCApKaJobSettings.reference_pair_molecules
        )

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
        assert len(captured["submissions"]) == 1
        # In "sub ... pka batch", jobs are not executed; only submission scripts are
        # generated, so reference molecules are not built at this stage.
        assert reference_pair_call_count["count"] == 0

    def test_sub_orca_pka_batch_first_exchange_rest_direct(
        self, tmp_path, monkeypatch
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule

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

        _attach_fake_array_server(monkeypatch, tmp_path, captured)
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
        assert len(captured["submissions"]) == 1
        children = _submitted_pka_child_jobs(captured)
        cli_args = _submitted_cli_args(captured)

        assert children[0].settings.scheme == "proton exchange"
        assert children[1].settings.scheme == "direct"
        assert len(cli_args) == 2
        assert "submit" in cli_args[0]
        assert "submit" in cli_args[1]
        assert "--reference" in cli_args[0]
        assert "--reference" not in cli_args[1]
        assert "--scheme" in cli_args[1]
        scheme_idx = cli_args[1].index("--scheme")
        assert cli_args[1][scheme_idx + 1] == "direct"

    def test_run_gaussian_pka_help_is_submission_only(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        _require_backend_pka_subcommand(run, "gaussian")
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
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
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        _require_backend_pka_subcommand(run, "orca")
        config_root = _write_test_backend_project(tmp_path, "orca")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
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

    def test_run_pka_help_keeps_output_analysis_commands(self):
        runner = CliRunner()
        result = runner.invoke(
            run, ["--no-scratch", "--fake", "pka", "--help"]
        )

        assert result.exit_code == 0, result.output
        assert "analyze" in result.output
        assert "batch-analyze" in result.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_auto_routes_to_batch(
        self, tmp_path, monkeypatch, backend
    ):
        """Table -f input should use batch workflow without requiring -pi."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
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
        assert "proton-index is required" not in result.output
        assert len(captured["submissions"]) == 1
        assert len(_submitted_pka_child_jobs(captured)) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_without_batch_subcommand(
        self, tmp_path, monkeypatch, backend
    ):
        """Omitting the batch subcommand still routes table input to batch."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "proton-index is required" not in result.output
        assert len(captured["submissions"]) == 1
        assert len(_submitted_pka_child_jobs(captured)) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_submit_subcommand_routes_to_batch(
        self, tmp_path, monkeypatch, backend
    ):
        """Explicit submit with table -f still uses row-wise batch processing."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "submit",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "proton-index is required" not in result.output
        assert len(captured["submissions"]) == 1
        assert len(_submitted_pka_child_jobs(captured)) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_batch_reconstructed_run_args_accept_proton_index(
        self, tmp_path, monkeypatch, backend
    ):
        """Batch run script should replay the full pka batch command under run."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
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
        assert captured["submissions"]

        cli_args = _submitted_cli_args(captured)
        assert len(cli_args) == 2
        assert all("submit" in args for args in cli_args)
        children = _submitted_pka_child_jobs(captured)
        assert len(children) == 2
        assert children[0].batch_entry["filepath"] in cli_args[0]
        assert children[1].batch_entry["filepath"] in cli_args[1]

        from chemsmart.jobs.gaussian.pka import GaussianpKaJob
        from chemsmart.jobs.orca.pka import ORCApKaJob

        run_labels = []

        def _fake_pka_run(self, **kwargs):
            run_labels.append(self.label)

        monkeypatch.setattr(GaussianpKaJob, "run", _fake_pka_run)
        monkeypatch.setattr(ORCApKaJob, "run", _fake_pka_run)

        for task_cli in cli_args:
            run_result = runner.invoke(
                run, ["--no-scratch", "--fake"] + list(task_cli)
            )
            assert run_result.exit_code == 0, run_result.output
            assert "proton-index is required" not in run_result.output
        assert len(run_labels) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_uses_coloured_proton_fragments(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """CDXML batch should create one job per coloured-proton fragment."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"submissions": []}
        _attach_fake_array_server(monkeypatch, tmp_path, captured)

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "Expected 5 fields" not in result.output
        assert "proton-index is required" not in result.output
        labels = [job.label for job in _submitted_pka_child_jobs(captured)]
        assert len(labels) == 2
        assert all("_frag" in label for label in labels)

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_uses_molecule_charge_without_parent_flags(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """CDXML batch should infer charge/mult from parsed Molecule objects."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"submissions": []}
        _attach_fake_array_server(monkeypatch, tmp_path, captured)

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        child_jobs = _submitted_pka_child_jobs(captured)
        assert len(child_jobs) == 2
        for job in child_jobs:
            assert job.settings.charge == 0
            assert job.settings.multiplicity == 1

    def test_get_pka_molecules_auto_assigns_charge_and_multiplicity(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.io.file import PKaCDXFile

        pka_mol = PKaCDXFile(colored_proton_cdxml_file).get_pka_molecules(
            index="-1"
        )
        assert pka_mol.charge == 0
        assert pka_mol.multiplicity == 1

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_reconstructed_scripts_target_single_fragment(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """Batch run script should replay the full CDXML batch command."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"submissions": []}
        _attach_fake_array_server(monkeypatch, tmp_path, captured)

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 1

        cli_args = _submitted_cli_args(captured)
        assert len(cli_args) == 2
        assert all("submit" in args for args in cli_args)
        assert all(
            str(colored_proton_two_molecule_cdxml_file) in args
            for args in cli_args
        )
        assert "--index" in cli_args[0]
        assert cli_args[0][cli_args[0].index("--index") + 1] == "1"
        assert "--index" in cli_args[1]
        assert cli_args[1][cli_args[1].index("--index") + 1] == "2"

        from chemsmart.jobs.gaussian.pka import GaussianpKaJob
        from chemsmart.jobs.orca.pka import ORCApKaJob

        run_labels = []

        def _fake_pka_run(self, **kwargs):
            run_labels.append(self.label)

        monkeypatch.setattr(GaussianpKaJob, "run", _fake_pka_run)
        monkeypatch.setattr(ORCApKaJob, "run", _fake_pka_run)

        for task_cli in cli_args:
            run_result = runner.invoke(
                run, ["--no-scratch", "--fake"] + list(task_cli)
            )
            assert run_result.exit_code == 0, run_result.output
            assert "proton-index is required" not in run_result.output
        assert len(run_labels) == 2
        assert all("_frag" in label for label in run_labels)

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_ignores_sibling_csv(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """CDXML batch must not fall back to a sibling CSV submission table."""
        _require_backend_pka_subcommand(sub, backend)
        sibling_csv = Path(colored_proton_two_molecule_cdxml_file).with_suffix(
            ".csv"
        )
        sibling_csv.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            "only_one_row.xyz,1,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"submissions": []}
        _attach_fake_array_server(monkeypatch, tmp_path, captured)

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        labels = [job.label for job in _submitted_pka_child_jobs(captured)]
        assert len(labels) == 2
        assert all("_frag" in label for label in labels)

    def test_pka_resolve_proton_index_accepts_explicit_index(self):
        from chemsmart.cli.pka import resolve_proton_index

        proton_index, molecules = resolve_proton_index("acid.xyz", 8, None)
        assert proton_index == 8
        assert molecules is None

    def test_resolve_pka_batch_row_auto_detects_coloured_proton(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row
        from chemsmart.io.molecules.structure import PKaMolecule

        proton_index, molecule = resolve_pka_batch_row(
            colored_proton_cdxml_file, proton_index=None
        )
        assert proton_index == 8
        assert isinstance(molecule, PKaMolecule)
        assert molecule.proton_index == 8

    def test_resolve_pka_batch_row_explicit_index_overrides_cdxml(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row
        from chemsmart.io.molecules.structure import Molecule

        proton_index, molecule = resolve_pka_batch_row(
            colored_proton_cdxml_file, proton_index=8
        )
        assert proton_index == 8
        assert isinstance(molecule, Molecule)

    def test_resolve_pka_batch_row_rejects_multi_molecule_cdxml(
        self, colored_proton_two_molecule_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row

        with pytest.raises(ValueError, match="single-molecule CDXML"):
            resolve_pka_batch_row(
                colored_proton_two_molecule_cdxml_file, proton_index=None
            )

    def test_resolve_pka_batch_row_requires_proton_index_for_xyz(
        self, single_molecule_xyz_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row

        with pytest.raises(ValueError, match="Missing proton_index"):
            resolve_pka_batch_row(single_molecule_xyz_file, proton_index=None)

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_cdxml_blank_proton_index_auto_detects(
        self, tmp_path, monkeypatch, backend, colored_proton_cdxml_file
    ):
        """CDXML rows with blank proton_index auto-detect the coloured proton."""
        _require_backend_pka_subcommand(sub, backend)

        table = tmp_path / "pka_cdxml.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{colored_proton_cdxml_file},,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"submissions": []}
        _attach_fake_array_server(monkeypatch, tmp_path, captured)

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
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
        assert len(captured["submissions"]) == 1
        children = _submitted_pka_child_jobs(captured)
        assert len(children) == 1
        assert children[0].settings.proton_index == 8

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_cdxml_explicit_proton_index_overrides(
        self, tmp_path, monkeypatch, backend, colored_proton_cdxml_file
    ):
        """Explicit table proton_index overrides CDXML coloured-proton detection."""
        _require_backend_pka_subcommand(sub, backend)

        table = tmp_path / "pka_cdxml.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{colored_proton_cdxml_file},8,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"submissions": []}
        _attach_fake_array_server(monkeypatch, tmp_path, captured)

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
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
        assert len(captured["submissions"]) == 1
        children = _submitted_pka_child_jobs(captured)
        assert len(children) == 1
        assert children[0].settings.proton_index == 8

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_rejects_multi_molecule_cdxml_row(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """Multi-molecule CDXML paths in a table row must fail clearly."""
        _require_backend_pka_subcommand(sub, backend)

        table = tmp_path / "pka_cdxml.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{colored_proton_two_molecule_cdxml_file},,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        fake_server.submit = lambda job, test=False, cli_args=None, **kw: None
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
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

        assert result.exit_code != 0
        assert "single-molecule CDXML" in result.output

    def test_orca_pka_job_generates_ha_and_a_subjobs(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch
    ):
        """ORCA pKa should prepare HA/A opt and SP jobs with Gaussian-style labels."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.opt import ORCAOptJob
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings
        from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )

        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )

        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="1a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )

        assert len(job.opt_jobs) == 2
        assert isinstance(job.protonated_job, ORCAOptJob)
        assert isinstance(job.conjugate_base_job, ORCAOptJob)
        assert job.protonated_job.label == "1a_pka_HA_opt"
        assert job.conjugate_base_job.label == "1a_pka_A_opt"
        assert job.conjugate_base_job.settings.charge == -1

        assert len(job.sp_jobs) == 2
        assert isinstance(job.protonated_sp_job, ORCASinglePointJob)
        assert isinstance(job.conjugate_base_sp_job, ORCASinglePointJob)
        assert job.protonated_sp_job.label == "1a_pka_HA_sp"
        assert job.conjugate_base_sp_job.label == "1a_pka_A_sp"

    def test_orca_pka_subjob_is_complete_uses_parent_folder(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch, tmp_path
    ):
        """Sub-jobs should detect completed outputs in the parent pKa folder."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )
        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="5a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)

        for name in ("5a_pka_HA_opt", "5a_pka_A_opt"):
            (tmp_path / f"{name}.out").write_text(
                "****ORCA TERMINATED NORMALLY****\n"
            )

        assert all(j.is_complete() for j in job.opt_jobs)

    def test_orca_pka_run_sp_jobs_after_completed_opt(
        self,
        single_molecule_xyz_file,
        orca_jobrunner_no_scratch,
        tmp_path,
        monkeypatch,
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )
        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="5a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)

        for name in ("5a_pka_HA_opt", "5a_pka_A_opt"):
            (tmp_path / f"{name}.out").write_text(
                "****ORCA TERMINATED NORMALLY****\n"
            )

        captured = {"sp_labels": []}

        def _fake_run_phase_jobs(*, jobs=None, jobs_factory=None, **kwargs):
            phase_jobs = jobs_factory() if jobs_factory is not None else jobs
            for child_job in phase_jobs:
                captured["sp_labels"].append(child_job.label)

        monkeypatch.setattr(
            "chemsmart.jobs.orca.pka.run_phase_jobs",
            _fake_run_phase_jobs,
        )
        monkeypatch.setattr(job, "_run_opt_jobs", lambda: None)
        monkeypatch.setattr(
            job, "_subjob_output", lambda *args, **kwargs: None
        )

        job._run()
        assert all(j.is_complete() for j in job.opt_jobs)
        assert captured["sp_labels"] == ["5a_pka_HA_sp", "5a_pka_A_sp"]

    def test_orca_pka_subjob_is_complete_recognizes_legacy_output(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch, tmp_path
    ):
        """Pre-rename ORCA pKa outputs should still count as complete."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )
        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="1a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)

        legacy_out = tmp_path / "1a_pka.out"
        legacy_out.write_text("****ORCA TERMINATED NORMALLY****\n")

        assert job._subjob_is_complete(
            job.protonated_job, legacy_label="1a_pka"
        )

    def test_orca_pka_run_executes_ha_and_a_opt_jobs(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch, monkeypatch
    ):
        """ORCA pKa opt phase should run both acid and conjugate-base jobs."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )

        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="1a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )

        captured = {"labels": []}

        def _fake_run_phase_jobs(*, jobs, **kwargs):
            for child_job in jobs:
                captured["labels"].append(child_job.label)

        monkeypatch.setattr(
            "chemsmart.jobs.orca.pka.run_phase_jobs",
            _fake_run_phase_jobs,
        )

        job._run_opt_jobs()
        assert captured["labels"] == ["1a_pka_HA_opt", "1a_pka_A_opt"]

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_run_pka_batch_table_processing(
        self, tmp_path, monkeypatch, backend
    ):
        """pKa table batch returns a BatchJob; run executes all child jobs."""
        _require_backend_pka_subcommand(run, backend)
        table = _build_pka_batch_table(tmp_path)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured = {"runs": []}

        from chemsmart.jobs.gaussian.batch import GaussianBatchJob
        from chemsmart.jobs.orca.batch import OrcaBatchJob

        def _fake_batch_run(self, **kwargs):
            for child in self.jobs:
                captured["runs"].append(child.label)

        monkeypatch.setattr(GaussianBatchJob, "run", _fake_batch_run)
        monkeypatch.setattr(OrcaBatchJob, "run", _fake_batch_run)

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
    def test_run_pka_batch_with_no_scratch(
        self, tmp_path, monkeypatch, backend
    ):
        """Explicit --no-scratch should not require a scratch directory."""
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

        from chemsmart.jobs.gaussian.batch import GaussianBatchJob
        from chemsmart.jobs.orca.batch import OrcaBatchJob

        def _fake_batch_run(self, **kwargs):
            for child in self.jobs:
                captured["runs"].append(child.label)

        monkeypatch.setattr(GaussianBatchJob, "run", _fake_batch_run)
        monkeypatch.setattr(OrcaBatchJob, "run", _fake_batch_run)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--fake",
                "--no-scratch",
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

    def test_run_rejects_non_job_batch_payload(self, pbs_server):
        """Scheduler-style batch payloads that are not Job lists stay blocked."""

        from chemsmart.cli.run import process_pipeline
        from chemsmart.jobs.runner import JobRunner

        ctx = click.Context(run)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = JobRunner(server=pbs_server, fake=True)

        with pytest.raises(
            ValueError, match="Expected a list of Job instances"
        ):
            process_pipeline.__wrapped__(ctx, ["not-a-job", "also-not-a-job"])
