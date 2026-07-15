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
    fake_server.submit = lambda job, test=False, cli_args=None, **kw: captured[
        "submissions"
    ].append((job, test, cli_args))
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    return table, captured


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
        self, tmp_path, monkeypatch
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
        assert first_args.index("submit") < first_args.index("--proton-index")

    def test_sub_orca_pka_batch_shared_reference_loaded_once(
        self, tmp_path, monkeypatch
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

    def test_sub_orca_pka_batch_first_exchange_rest_direct(
        self, tmp_path, monkeypatch
    ):
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
        assert (
            first_args[first_args.index("--scheme") + 1] == "proton exchange"
        )
        assert "--reference" in first_args

        assert "--scheme" in second_args
        assert second_args[second_args.index("--scheme") + 1] == "direct"
        assert "--reference" not in second_args
        assert "--reference-proton-index" not in second_args
        assert "--reference-charge" not in second_args
        assert "--reference-multiplicity" not in second_args

    def test_run_gaussian_pka_help_is_submission_only(
        self, single_molecule_xyz_file
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
        self, single_molecule_xyz_file
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
        assert len(captured["submissions"]) == 2

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
        assert len(captured["submissions"]) == 2

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
        assert len(captured["submissions"]) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_batch_reconstructed_run_args_accept_proton_index(
        self, tmp_path, monkeypatch, backend
    ):
        """Per-row chemsmart_run_*.py args must parse --proton-index under run."""
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

        cli_args = captured["submissions"][0][2]
        assert "submit" in cli_args
        assert "--proton-index" in cli_args
        assert cli_args.index("submit") < cli_args.index("--proton-index")

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            return None

        monkeypatch.setattr(Job, "run", _fake_run)

        run_result = runner.invoke(run, ["--no-scratch", "--fake"] + cli_args)
        assert run_result.exit_code == 0, run_result.output
        assert "proton-index is required" not in run_result.output

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

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured = {"labels": []}
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "labels"
            ].append(job.label)
        )
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
        assert len(captured["labels"]) == 2
        assert all("_frag" in label for label in captured["labels"])

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

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured = {"jobs": []}
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "jobs"
            ].append(job)
        )
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
                colored_proton_two_molecule_cdxml_file,
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["jobs"]) == 2
        for job in captured["jobs"]:
            assert job.settings.charge == 0
            assert job.settings.multiplicity == 1
            assert job._batch_entry["charge"] == 0
            assert job._batch_entry["multiplicity"] == 1
            assert job._batch_entry["label"] == job.label

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
        """Each CDXML fragment script must submit only that fragment, not re-batch all."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured = {"submissions": []}
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
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
        assert len(captured["submissions"]) == 2

        fragment_indices = []
        for job, _test, cli_args in captured["submissions"]:
            assert "batch" not in cli_args
            assert "submit" in cli_args
            assert "--proton-index" in cli_args
            assert "--index" in cli_args
            assert "--label" in cli_args
            assert cli_args[cli_args.index("--label") + 1] == job.label
            fragment_indices.append(cli_args[cli_args.index("--index") + 1])

        assert fragment_indices == ["1", "2"]

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            return None

        monkeypatch.setattr(Job, "run", _fake_run)

        for job, _test, cli_args in captured["submissions"]:
            run_labels = []

            def _fake_run(self):
                run_labels.append(self.label)
                return None

            monkeypatch.setattr(Job, "run", _fake_run)
            run_result = runner.invoke(
                run, ["--no-scratch", "--fake"] + cli_args
            )
            assert run_result.exit_code == 0, run_result.output
            assert "proton-index is required" not in run_result.output
            assert run_labels == [job.label]

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

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured = {"labels": []}
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "labels"
            ].append(job.label)
        )
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
        assert len(captured["labels"]) == 2
        assert all("_frag" in label for label in captured["labels"])

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

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured = {"submissions": []}
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
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

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 1
        job = captured["submissions"][0][0]
        assert job.settings.proton_index == 8

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

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured = {"submissions": []}
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
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

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 1
        job = captured["submissions"][0][0]
        assert job.settings.proton_index == 8

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
    def test_run_pka_batch_defaults_to_no_scratch(
        self, tmp_path, monkeypatch, backend
    ):
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

    def test_run_rejects_non_job_batch_payload(self, pbs_server):
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
