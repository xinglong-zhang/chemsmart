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

    def test_compute_pka_direct_requires_delta_g_proton(self):
        from chemsmart.cli.pka import compute_pka

        with pytest.raises(ValueError, match="delta_G_proton is required"):
            compute_pka(
                ha_gas_file="ha.log",
                a_gas_file="a.log",
                ha_solv_file="has.log",
                a_solv_file="as.log",
                scheme="direct",
                delta_G_proton=None,
            )

    def test_compute_pka_proton_exchange_requires_pka_reference(self):
        from chemsmart.cli.pka import compute_pka

        with pytest.raises(ValueError, match="pka_reference is required"):
            compute_pka(
                ha_gas_file="ha.log",
                a_gas_file="a.log",
                scheme="proton exchange",
                pka_reference=None,
            )

    def test_compute_pka_proton_exchange_requires_all_files(self):
        from chemsmart.cli.pka import compute_pka

        with pytest.raises(ValueError, match="Missing required files"):
            compute_pka(
                ha_gas_file="ha.log",
                a_gas_file="a.log",
                scheme="proton exchange",
                pka_reference=6.75,
                # href_gas_file, ref_gas_file, ha_solv_file, a_solv_file,
                # href_solv_file, ref_solv_file all left as None
            )

    def test_compute_pka_direct_requires_solvent_files(self):
        from chemsmart.cli.pka import compute_pka

        with pytest.raises(
            ValueError, match="ha_solv_file and a_solv_file are required"
        ):
            compute_pka(
                ha_gas_file="ha.log",
                a_gas_file="a.log",
                scheme="direct",
                delta_G_proton=-265.9,
                # ha_solv_file/a_solv_file omitted
            )

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

    def test_print_pka_summary_direct_scheme(
        self,
        capsys,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        gaussian_pKa_HA_single_point_outputfile,
        gaussian_pKa_A_single_point_outputfile,
    ):
        """The 'direct' scheme branch of print_pka_summary is otherwise
        only reached via print_pka_summary mocked out in CLI tests, so
        exercise the real formatted-output body directly here."""
        from chemsmart.cli.pka import print_pka_summary

        print_pka_summary(
            ha_gas_file=gaussian_pKa_HA_optimization_outputfile,
            a_gas_file=gaussian_pKa_A_optimization_outputfile,
            ha_solv_file=gaussian_pKa_HA_single_point_outputfile,
            a_solv_file=gaussian_pKa_A_single_point_outputfile,
            scheme="direct",
            delta_G_proton=-265.9,
            temperature=373.15,
        )
        output = capsys.readouterr().out
        assert "Direct Dissociation Scheme" in output
        assert "Computed pKa(HA)" in output

    def test_compute_pka_thermochemistry_includes_href_and_ref(
        self,
        gaussian_pKa_HA_optimization_outputfile,
        gaussian_pKa_A_optimization_outputfile,
        gaussian_pKa_HB_optimization_outputfile,
        gaussian_pKa_B_optimization_outputfile,
    ):
        """HRef/Ref branches (only reached with a 4-file proton-exchange
        scheme) are not exercised by the HA/A-only tests elsewhere."""
        from chemsmart.cli.pka import compute_pka_thermochemistry

        results = compute_pka_thermochemistry(
            ha_file=gaussian_pKa_HA_optimization_outputfile,
            a_file=gaussian_pKa_A_optimization_outputfile,
            href_file=gaussian_pKa_HB_optimization_outputfile,
            ref_file=gaussian_pKa_B_optimization_outputfile,
        )
        assert results["HA"]["name"] == "HA"
        assert results["A"]["name"] == "A-"
        assert results["HRef"]["name"] == "HRef"
        assert results["Ref"]["name"] == "Ref-"

    def test_require_pka_charge_multiplicity_reports_missing_charge(self):
        from types import SimpleNamespace

        from chemsmart.cli.pka import require_pka_charge_multiplicity

        opt_settings = SimpleNamespace(charge=None, multiplicity=1)
        with pytest.raises(click.UsageError, match="-c/--charge"):
            require_pka_charge_multiplicity(opt_settings)

    def test_require_pka_charge_multiplicity_reports_missing_multiplicity(
        self,
    ):
        from types import SimpleNamespace

        from chemsmart.cli.pka import require_pka_charge_multiplicity

        opt_settings = SimpleNamespace(charge=0, multiplicity=None)
        with pytest.raises(click.UsageError, match="-m/--multiplicity"):
            require_pka_charge_multiplicity(opt_settings)

    def test_require_pka_charge_multiplicity_includes_source_hint(self):
        from types import SimpleNamespace

        from chemsmart.cli.pka import require_pka_charge_multiplicity

        opt_settings = SimpleNamespace(charge=None, multiplicity=None)
        with pytest.raises(click.UsageError, match=r"\(from row 2\)"):
            require_pka_charge_multiplicity(
                opt_settings, source_hint="from row 2"
            )

    def test_require_pka_charge_multiplicity_passes_when_both_set(self):
        from types import SimpleNamespace

        from chemsmart.cli.pka import require_pka_charge_multiplicity

        opt_settings = SimpleNamespace(charge=0, multiplicity=1)
        assert require_pka_charge_multiplicity(opt_settings) is None

    def test_is_pka_batch_invocation_false_for_non_pka_subcommand(self):
        from types import SimpleNamespace

        from chemsmart.cli.pka import is_pka_batch_invocation

        ctx = SimpleNamespace(invoked_subcommand="other", parent=None)
        assert is_pka_batch_invocation(ctx) is False

    def test_is_pka_batch_invocation_false_when_submit_present(self):
        from types import SimpleNamespace

        from chemsmart.cli.pka import is_pka_batch_invocation

        ctx = SimpleNamespace(
            invoked_subcommand="pka",
            args=["batch", "submit"],
            parent=None,
        )
        assert is_pka_batch_invocation(ctx) is False

    def test_is_pka_batch_invocation_true_when_batch_present(self):
        from types import SimpleNamespace

        from chemsmart.cli.pka import is_pka_batch_invocation

        ctx = SimpleNamespace(
            invoked_subcommand="pka",
            args=["batch"],
            parent=None,
        )
        assert is_pka_batch_invocation(ctx) is True

    def test_compute_pka_thermochemistry_href_ref_only(
        self,
        gaussian_pKa_HB_optimization_outputfile,
        gaussian_pKa_B_optimization_outputfile,
    ):
        """ha_file/a_file may be omitted when only reference species are
        needed, exercising the "not provided" branches for each."""
        from chemsmart.cli.pka import compute_pka_thermochemistry

        results = compute_pka_thermochemistry(
            href_file=gaussian_pKa_HB_optimization_outputfile,
            ref_file=gaussian_pKa_B_optimization_outputfile,
        )
        assert "HA" not in results
        assert "A" not in results
        assert results["HRef"]["name"] == "HRef"
        assert results["Ref"]["name"] == "Ref-"

    def test_resolve_pka_analysis_scheme_direct_requires_delta_g(self):
        from chemsmart.cli.pka import _resolve_pka_analysis_scheme

        with pytest.raises(
            click.UsageError, match="delta-g-proton is required"
        ):
            _resolve_pka_analysis_scheme(scheme="direct", delta_g_proton=None)

    def test_resolve_pka_analysis_scheme_ignores_delta_g_for_exchange(
        self, caplog
    ):
        from chemsmart.cli.pka import _resolve_pka_analysis_scheme

        with caplog.at_level("INFO"):
            scheme = _resolve_pka_analysis_scheme(
                scheme="proton exchange", delta_g_proton=-265.9
            )
        assert scheme == "proton exchange"
        assert "Ignoring -dG/--delta-g-proton" in caplog.text

    def test_resolve_pka_analysis_scheme_defaults_to_proton_exchange(self):
        from chemsmart.cli.pka import _resolve_pka_analysis_scheme

        assert (
            _resolve_pka_analysis_scheme(scheme=None, delta_g_proton=None)
            == "proton exchange"
        )

    def test_validate_direct_analyze_files_reports_missing(self, tmp_path):
        from chemsmart.cli.pka import validate_direct_analyze_files

        with pytest.raises(click.UsageError, match="all four output files"):
            validate_direct_analyze_files(
                ha=str(tmp_path / "ha.log"), a=None, ha_solv=None, a_solv=None
            )

    def test_validate_direct_analyze_files_reports_not_found(self, tmp_path):
        from chemsmart.cli.pka import validate_direct_analyze_files

        ha = tmp_path / "ha.log"
        ha.write_text("Gaussian, Inc.\n")
        missing_a = tmp_path / "missing_a.log"
        with pytest.raises(click.UsageError, match="File not found"):
            validate_direct_analyze_files(
                ha=str(ha),
                a=str(missing_a),
                ha_solv=str(ha),
                a_solv=str(ha),
            )

    def test_validate_direct_analyze_files_passes_when_all_present(
        self, tmp_path
    ):
        from chemsmart.cli.pka import validate_direct_analyze_files

        ha = tmp_path / "ha.log"
        ha.write_text("Gaussian, Inc.\n")
        assert (
            validate_direct_analyze_files(
                ha=str(ha), a=str(ha), ha_solv=str(ha), a_solv=str(ha)
            )
            is None
        )

    def test_auto_discover_direct_pka_files_reports_missing_companions(
        self, tmp_path
    ):
        from chemsmart.cli.pka import _auto_discover_direct_pka_files

        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        with pytest.raises(click.UsageError, match="Auto-discovery"):
            _auto_discover_direct_pka_files(str(ha_gas))

    def test_auto_discover_direct_pka_files_finds_all_companions(
        self, tmp_path
    ):
        """Success path (program auto-detected from the HA gas file's
        signature) when every companion output already exists on disk."""
        from chemsmart.cli.pka import _auto_discover_direct_pka_files

        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"target{suffix}.log").write_text("Gaussian, Inc.\n")

        results = _auto_discover_direct_pka_files(str(ha_gas))
        assert all(key in results for key in ("a", "ha_solv", "a_solv"))

    def test_validate_reference_options_noop_when_no_reference(self):
        from chemsmart.cli.pka import validate_reference_options

        shared = {"reference": None, "scheme": "direct"}
        assert validate_reference_options(shared) is None

    def test_validate_reference_options_requires_proton_exchange_scheme(self):
        from chemsmart.cli.pka import validate_reference_options

        shared = {"reference": "ref.xyz", "scheme": "direct"}
        with pytest.raises(click.UsageError, match="proton exchange"):
            validate_reference_options(shared)

    def test_validate_reference_options_reports_missing_fields(self):
        from chemsmart.cli.pka import validate_reference_options

        shared = {
            "reference": "ref.xyz",
            "scheme": "proton exchange",
            "reference_proton_index": None,
            "reference_color_code": None,
            "reference_charge": None,
            "reference_multiplicity": None,
        }
        with pytest.raises(click.UsageError, match="reference-proton-index"):
            validate_reference_options(shared)

    def test_validate_reference_options_passes_when_all_set(self):
        from chemsmart.cli.pka import validate_reference_options

        shared = {
            "reference": "ref.xyz",
            "scheme": "proton exchange",
            "reference_proton_index": 8,
            "reference_color_code": None,
            "reference_charge": 0,
            "reference_multiplicity": 1,
        }
        assert validate_reference_options(shared) is None

    def test_is_existing_output_path_variants(self, tmp_path):
        from chemsmart.cli.pka import _is_existing_output_path

        assert _is_existing_output_path(None) is False
        assert _is_existing_output_path("") is False
        real_file = tmp_path / "a.log"
        real_file.write_text("x")
        assert _is_existing_output_path(str(real_file)) is True
        assert _is_existing_output_path(str(tmp_path / "missing.log")) is False

    def test_validate_pka_table_program_flags_mismatch(self, tmp_path):
        from chemsmart.cli.pka import _validate_pka_table_program

        orca_file = tmp_path / "a.out"
        orca_file.write_text("* O   R   C   A *\n")

        class _FakeTable:
            entries = [{"ha_gas": str(orca_file)}]

        with pytest.raises(click.UsageError, match="was detected as"):
            _validate_pka_table_program(_FakeTable(), "gaussian")

    def test_validate_pka_table_program_passes_when_matching(self, tmp_path):
        from chemsmart.cli.pka import _validate_pka_table_program

        orca_file = tmp_path / "a.out"
        orca_file.write_text("* O   R   C   A *\n")

        class _FakeTable:
            entries = [{"ha_gas": str(orca_file)}]

        assert _validate_pka_table_program(_FakeTable(), "orca") is None

    def test_validate_analyze_files_reports_missing_gas_and_solv(self):
        from chemsmart.cli.pka import validate_analyze_files

        with pytest.raises(click.UsageError, match="all 8 output files"):
            validate_analyze_files(
                ha=None,
                a=None,
                href=None,
                ref=None,
                ha_solv=None,
                a_solv=None,
                href_solv=None,
                ref_solv=None,
                reference_pka=6.75,
            )

    def test_validate_analyze_files_requires_reference_pka(self, tmp_path):
        from chemsmart.cli.pka import validate_analyze_files

        f = tmp_path / "x.log"
        f.write_text("Gaussian, Inc.\n")
        with pytest.raises(click.UsageError, match="reference-pka"):
            validate_analyze_files(
                ha=str(f),
                a=str(f),
                href=str(f),
                ref=str(f),
                ha_solv=str(f),
                a_solv=str(f),
                href_solv=str(f),
                ref_solv=str(f),
                reference_pka=None,
            )

    def test_validate_analyze_files_reports_nonexistent_paths(self, tmp_path):
        from chemsmart.cli.pka import validate_analyze_files

        f = tmp_path / "x.log"
        f.write_text("Gaussian, Inc.\n")
        missing = str(tmp_path / "missing.log")
        with pytest.raises(click.UsageError, match="do not exist"):
            validate_analyze_files(
                ha=missing,
                a=str(f),
                href=str(f),
                ref=str(f),
                ha_solv=str(f),
                a_solv=str(f),
                href_solv=str(f),
                ref_solv=str(f),
                reference_pka=6.75,
            )

    def test_validate_analyze_files_reports_nonexistent_solv_path(
        self, tmp_path
    ):
        from chemsmart.cli.pka import validate_analyze_files

        f = tmp_path / "x.log"
        f.write_text("Gaussian, Inc.\n")
        missing_solv = str(tmp_path / "missing_solv.log")
        with pytest.raises(click.UsageError, match="do not exist"):
            validate_analyze_files(
                ha=str(f),
                a=str(f),
                href=str(f),
                ref=str(f),
                ha_solv=missing_solv,
                a_solv=str(f),
                href_solv=str(f),
                ref_solv=str(f),
                reference_pka=6.75,
            )

    def test_validate_analyze_files_passes_when_all_present(self, tmp_path):
        from chemsmart.cli.pka import validate_analyze_files

        f = tmp_path / "x.log"
        f.write_text("Gaussian, Inc.\n")
        assert (
            validate_analyze_files(
                ha=str(f),
                a=str(f),
                href=str(f),
                ref=str(f),
                ha_solv=str(f),
                a_solv=str(f),
                href_solv=str(f),
                ref_solv=str(f),
                reference_pka=6.75,
            )
            is None
        )

    def test_auto_discover_pka_files_reports_missing_companions(
        self, tmp_path
    ):
        from chemsmart.cli.pka import _auto_discover_pka_files

        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        href_gas = tmp_path / "ref_HRef_opt.log"
        href_gas.write_text("Gaussian, Inc.\n")
        with pytest.raises(click.UsageError, match="Auto-discovery"):
            _auto_discover_pka_files(
                str(ha_gas), str(href_gas), program="gaussian"
            )

    def test_auto_discover_pka_files_finds_all_companions_auto_program(
        self, tmp_path
    ):
        """Success path with program=None (auto-detected via
        get_program_type_from_file from the HA gas file's signature)."""
        from chemsmart.cli.pka import _auto_discover_pka_files

        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"target{suffix}.log").write_text("Gaussian, Inc.\n")

        href_gas = tmp_path / "ref_HRef_opt.log"
        href_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_Ref_opt", "_HRef_sp", "_Ref_sp"):
            (tmp_path / f"ref{suffix}.log").write_text("Gaussian, Inc.\n")

        results = _auto_discover_pka_files(str(ha_gas), str(href_gas))
        assert all(
            key in results
            for key in (
                "a",
                "ha_solv",
                "a_solv",
                "ref",
                "href_solv",
                "ref_solv",
            )
        )

    def test_run_pka_analyze_direct_auto_discovers_companion_files(
        self, tmp_path, monkeypatch
    ):
        """`-s direct analyze -ha ...` alone should auto-discover
        a/ha-solv/a-solv from the standard suffix convention."""
        monkeypatch.chdir(tmp_path)
        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"target{suffix}.log").write_text("Gaussian, Inc.\n")

        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "-s",
                "direct",
                "-dG",
                "-265.9",
                "analyze",
                "-ha",
                str(ha_gas),
            ],
        )

        assert result.exit_code == 0, result.output
        assert called["kwargs"]["ha_gas_file"] == str(ha_gas)
        assert called["kwargs"]["a_gas_file"] == str(
            tmp_path / "target_pka_A_opt.log"
        )
        assert called["kwargs"]["ha_solv_file"] == str(
            tmp_path / "target_pka_HA_sp.log"
        )

    def test_run_pka_analyze_proton_exchange_auto_discovers_companion_files(
        self, tmp_path, monkeypatch
    ):
        """`analyze -ha ... -hr ...` (default scheme) should auto-discover
        the remaining six companion output files."""
        monkeypatch.chdir(tmp_path)
        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"target{suffix}.log").write_text("Gaussian, Inc.\n")

        href_gas = tmp_path / "ref_HRef_opt.log"
        href_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_Ref_opt", "_HRef_sp", "_Ref_sp"):
            (tmp_path / f"ref{suffix}.log").write_text("Gaussian, Inc.\n")

        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "analyze",
                "-ha",
                str(ha_gas),
                "-hr",
                str(href_gas),
                "-rp",
                "6.75",
            ],
        )

        assert result.exit_code == 0, result.output
        assert called["kwargs"]["ha_gas_file"] == str(ha_gas)
        assert called["kwargs"]["href_gas_file"] == str(href_gas)
        assert called["kwargs"]["a_gas_file"] == str(
            tmp_path / "target_pka_A_opt.log"
        )
        assert called["kwargs"]["ref_gas_file"] == str(
            tmp_path / "ref_Ref_opt.log"
        )

    def test_run_pka_batch_analyze_bad_table_becomes_usage_error(
        self, tmp_path
    ):
        """A table whose referenced companion files don't exist on disk
        should surface prepare()'s ValueError as a clean click.UsageError,
        not an uncaught exception."""
        table = tmp_path / "table.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            "target,,,,,ref_HA_opt.log,ref_A_opt.log,ref_HA_sp.log,ref_A_sp.log,10.6\n"
        )
        runner = CliRunner()
        result = runner.invoke(
            run,
            ["pka", "batch-analyze", "-o", str(table)],
        )
        assert result.exit_code != 0
        assert "File not found" in result.output

    def test_run_pka_batch_analyze_explicit_program_validates_table(
        self, tmp_path, monkeypatch
    ):
        """Explicit -p should be checked against every table file's
        detected program (not just left to 'auto')."""
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

        _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "batch-analyze",
                "-o",
                str(table),
                "-p",
                "gaussian",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "pKa" in result.output

    def test_run_pka_analyze_direct_requires_ha(self):
        runner = CliRunner()
        result = runner.invoke(
            run,
            ["pka", "-s", "direct", "-dG", "-265.9", "analyze"],
        )
        assert result.exit_code != 0
        assert "-ha/--ha is required" in result.output

    def test_run_pka_analyze_proton_exchange_skips_discovery_without_ha_and_href(
        self,
    ):
        """Without both -ha and -hr, auto-discovery can't run at all, so
        this should fall straight through to the missing-files error."""
        runner = CliRunner()
        result = runner.invoke(run, ["pka", "analyze", "-rp", "6.75"])
        assert result.exit_code != 0
        assert "all 8 output files are required" in result.output

    def test_run_pka_analyze_direct_partial_auto_discovery(
        self, tmp_path, monkeypatch
    ):
        """Explicitly providing one optional file (-a) alongside -ha
        should still auto-discover only the remaining missing ones."""
        monkeypatch.chdir(tmp_path)
        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        explicit_a = tmp_path / "custom_a.log"
        explicit_a.write_text("Gaussian, Inc.\n")
        # _auto_discover_direct_pka_files unconditionally discovers and
        # checks all three companions (even ones already explicitly
        # provided), so the naming-convention "a" file must exist too.
        for suffix in ("_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"target{suffix}.log").write_text("Gaussian, Inc.\n")

        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "-s",
                "direct",
                "-dG",
                "-265.9",
                "analyze",
                "-ha",
                str(ha_gas),
                "-a",
                str(explicit_a),
            ],
        )
        assert result.exit_code == 0, result.output
        assert called["kwargs"]["a_gas_file"] == str(explicit_a)
        assert called["kwargs"]["ha_solv_file"] == str(
            tmp_path / "target_pka_HA_sp.log"
        )

    def test_run_pka_analyze_proton_exchange_partial_auto_discovery(
        self, tmp_path, monkeypatch
    ):
        """Explicitly providing one optional file (-a) alongside -ha/-hr
        should still auto-discover only the remaining missing ones."""
        monkeypatch.chdir(tmp_path)
        ha_gas = tmp_path / "target_pka_HA_opt.log"
        ha_gas.write_text("Gaussian, Inc.\n")
        explicit_a = tmp_path / "custom_a.log"
        explicit_a.write_text("Gaussian, Inc.\n")
        # _auto_discover_pka_files unconditionally discovers and checks
        # all six companions (even ones already explicitly provided), so
        # the naming-convention "a" file must exist too.
        for suffix in ("_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"target{suffix}.log").write_text("Gaussian, Inc.\n")

        href_gas = tmp_path / "ref_HRef_opt.log"
        href_gas.write_text("Gaussian, Inc.\n")
        for suffix in ("_Ref_opt", "_HRef_sp", "_Ref_sp"):
            (tmp_path / f"ref{suffix}.log").write_text("Gaussian, Inc.\n")

        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "analyze",
                "-ha",
                str(ha_gas),
                "-hr",
                str(href_gas),
                "-a",
                str(explicit_a),
                "-rp",
                "6.75",
            ],
        )
        assert result.exit_code == 0, result.output
        assert called["kwargs"]["a_gas_file"] == str(explicit_a)
        assert called["kwargs"]["ref_gas_file"] == str(
            tmp_path / "ref_Ref_opt.log"
        )

    def test_resolve_pka_submit_proton_options_color_code_from_parent(self):
        """When the parent group already resolved -cc, it must be used
        directly without falling back to ctx.obj."""
        from types import SimpleNamespace

        from chemsmart.cli.pka import resolve_pka_submit_proton_options

        parent = SimpleNamespace(
            params={"proton_index": None, "color_code": 4}
        )
        ctx = SimpleNamespace(
            parent=parent,
            obj={"pka_proton_index": None, "pka_color_code": None},
        )
        proton_index, color_code = resolve_pka_submit_proton_options(ctx)
        assert color_code == 4
        assert ctx.obj["pka_color_code"] == 4

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

    def test_run_gaussian_pka_no_subcommand_auto_dispatches_to_submit(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        """``pka`` with no explicit submit/batch subcommand (and a
        non-table input file) should auto-dispatch to submit()."""
        _require_backend_pka_subcommand(run, "gaussian")
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.jobs.job import Job

        monkeypatch.setattr(Job, "run", lambda self: None)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--no-scratch",
                "--fake",
                "gaussian",
                "-p",
                "test",
                "-c",
                "0",
                "-m",
                "1",
                "-f",
                single_molecule_xyz_file,
                "pka",
                "-s",
                "direct",
                "-dG",
                "-265.9",
                "-pi",
                "19",
            ],
        )
        assert result.exit_code == 0, result.output

    def test_gaussian_pka_batch_requires_filename(self, tmp_path, monkeypatch):
        """``pka batch`` without a parent -f/--filename (e.g. a
        --pubchem-only invocation) raises a clear UsageError."""
        from unittest.mock import MagicMock, patch

        _require_backend_pka_subcommand(run, "gaussian")
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[MagicMock()],
        ):
            runner = CliRunner()
            result = runner.invoke(
                run,
                [
                    "--no-scratch",
                    "--fake",
                    "gaussian",
                    "-p",
                    "test",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "-l",
                    "ammonia",
                    "--pubchem",
                    "222",
                    "pka",
                    "-s",
                    "direct",
                    "-dG",
                    "-265.9",
                    "batch",
                ],
            )
        assert result.exit_code != 0
        assert "Batch mode requires" in result.output

    def test_gaussian_pka_batch_malformed_table_becomes_usage_error(
        self, tmp_path, monkeypatch
    ):
        _require_backend_pka_subcommand(run, "gaussian")
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        table = tmp_path / "bad_table.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            "does_not_exist.xyz,2,0,1\n"
        )

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
                str(table),
                "pka",
                "-s",
                "direct",
                "-dG",
                "-265.9",
                "batch",
            ],
        )
        assert result.exit_code != 0

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

    def test_resolve_proton_index_color_code_requires_cdxml(
        self, single_molecule_xyz_file
    ):
        from chemsmart.cli.pka import resolve_proton_index

        with pytest.raises(ValueError, match="color-code can only be used"):
            resolve_proton_index(
                single_molecule_xyz_file, proton_index=None, color_code=2
            )

    def test_resolve_proton_index_detects_submission_table(self, tmp_path):
        from chemsmart.cli.pka import resolve_proton_index

        table = _build_pka_batch_table(tmp_path)
        with pytest.raises(ValueError, match="Use the 'batch' subcommand"):
            resolve_proton_index(str(table), proton_index=None)

    def test_resolve_proton_index_requires_index_for_plain_file(
        self, single_molecule_xyz_file
    ):
        from chemsmart.cli.pka import resolve_proton_index

        with pytest.raises(ValueError, match="-pi/--proton-index is required"):
            resolve_proton_index(single_molecule_xyz_file, proton_index=None)

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

    def test_resolve_pka_batch_row_wraps_color_detection_failure(
        self, colored_proton_cdxml_file
    ):
        """An unmatched -cc/--color-code should surface a clear,
        re-wrapped ValueError rather than the raw lookup failure."""
        from chemsmart.cli.pka import resolve_pka_batch_row

        with pytest.raises(
            ValueError, match="Could not auto-detect proton from CDXML colour"
        ):
            resolve_pka_batch_row(
                colored_proton_cdxml_file, proton_index=None, color_code=999
            )

    def test_batch_pka_jobs_from_cdxml_wraps_resolve_error_as_usage_error(
        self, single_molecule_xyz_file
    ):
        """resolve_proton_index's ValueError (e.g. -cc on a non-CDXML
        file) must surface as a click.UsageError, not an uncaught one."""
        from types import SimpleNamespace

        from chemsmart.cli.pka import batch_pka_jobs_from_cdxml

        ctx = SimpleNamespace(
            obj={
                "filename": single_molecule_xyz_file,
                "pka_shared": {},
                "pka_proton_index": None,
                "pka_color_code": 2,
            },
            parent=None,
        )
        with pytest.raises(click.UsageError, match="color-code can only"):
            batch_pka_jobs_from_cdxml(
                ctx,
                skip_completed=False,
                create_jobs_fn=lambda *a, **kw: None,
                invoke_submit_fn=lambda *a, **kw: None,
            )

    def test_batch_pka_jobs_from_cdxml_calls_invoke_submit_for_single_molecule(
        self, single_molecule_xyz_file
    ):
        """When resolve_proton_index resolves a single molecule (no
        per-fragment list), the submit path is invoked, not create_jobs."""
        from types import SimpleNamespace

        from chemsmart.cli.pka import batch_pka_jobs_from_cdxml

        ctx = SimpleNamespace(
            obj={
                "filename": single_molecule_xyz_file,
                "pka_shared": {},
                "pka_proton_index": 2,
                "pka_color_code": None,
            },
            parent=None,
        )
        captured = {}

        def _fake_invoke_submit(ctx, skip_completed, proton_index, color_code):
            captured["proton_index"] = proton_index
            return "submitted"

        result = batch_pka_jobs_from_cdxml(
            ctx,
            skip_completed=False,
            create_jobs_fn=lambda *a, **kw: pytest.fail(
                "create_jobs_fn should not be called for a single molecule"
            ),
            invoke_submit_fn=_fake_invoke_submit,
        )
        assert result == "submitted"
        assert captured["proton_index"] == 2

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

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            captured["runs"].append(self.label)

        monkeypatch.setattr(Job, "run", _fake_run)

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
            ValueError, match="Batch job submission is not supported"
        ):
            process_pipeline.__wrapped__(ctx, ["not-a-job", "also-not-a-job"])


def _build_pka_shared(**overrides):
    """Minimal 'pka_shared' dict matching the keys the `pka` group builds."""
    shared = dict(
        scheme="direct",
        reference=None,
        reference_proton_index=None,
        reference_color_code=None,
        reference_charge=None,
        reference_multiplicity=None,
        reference_conjugate_base_charge=None,
        reference_conjugate_base_multiplicity=None,
        delta_g_proton=None,
        conjugate_base_charge=None,
        conjugate_base_multiplicity=None,
        solvent_model=None,
        solvent_id=None,
        temperature=None,
        concentration=None,
        pressure=None,
        cutoff_entropy_grimme=None,
        cutoff_enthalpy=None,
        entropy_method=None,
        skip_completed=False,
    )
    shared.update(overrides)
    return shared


class TestPkaCliDirectBranchCoverage:
    """Targeted coverage for chemsmart/cli/gaussian/pka.py branches not
    reached by the higher-level batch/analyze tests above: submit()'s
    own multi-fragment CDXML routing, the multi-molecule-index job
    list, the proton-exchange reference validation in batch(), and the
    defensive job_settings-not-set branches in batch() and
    _create_pka_jobs_from_molecules()."""

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_submit_subcommand_handles_multifragment_cdxml(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """Explicit 'submit' (not 'batch') with a multi-fragment CDXML
        should still create one job per fragment, via submit()'s own
        `pka_molecules is not None` branch rather than
        batch_pka_jobs_from_cdxml."""
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
                "submit",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["labels"]) == 2
        assert all("_frag" in label for label in captured["labels"])

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_submit_multi_index_creates_one_job_per_index(
        self,
        tmp_path,
        monkeypatch,
        backend,
        two_rotated_molecules_xyz_file,
    ):
        """Selecting multiple molecule indices from a multi-structure
        file together with an explicit -pi proton index takes the
        `len(molecules) > 1 and molecule_indices` branch, creating one
        job per selected index instead of a single job."""
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
                two_rotated_molecules_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-i",
                "1,2",
                "pka",
                "-s",
                "direct",
                "-pi",
                "7",
                "submit",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["labels"]) == 2
        assert all("_idx" in label for label in captured["labels"])

    def test_batch_proton_exchange_missing_all_reference_options_raises(
        self, tmp_path, monkeypatch
    ):
        """scheme=proton exchange with no reference options at all should
        list every missing option in one UsageError."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        table = _build_pka_batch_table(tmp_path)

        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")
        batch_cmd = pka_mod.pka.commands["batch"]

        shared = _build_pka_shared(scheme="proton exchange")
        ctx = click.Context(batch_cmd)
        ctx.obj = {
            "pka_shared": shared,
            "filename": str(table),
            "jobrunner": None,
            "project_settings": project_settings,
            "job_settings": None,
            "keywords": {},
        }
        with ctx:
            with pytest.raises(click.UsageError) as exc_info:
                batch_cmd.callback(
                    skip_completed=False, proton_index=None, color_code=None
                )
        message = str(exc_info.value)
        assert "-r/--reference" in message
        assert "-rc/--reference-charge" in message
        assert "-rm/--reference-multiplicity" in message

    def test_batch_proton_exchange_noncdxml_reference_missing_proton_index_raises(
        self, tmp_path, monkeypatch
    ):
        """A non-CDXML reference file without --reference-proton-index
        cannot be auto-detected, so it must be reported as missing."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        table = _build_pka_batch_table(tmp_path)
        reference = tmp_path / "reference_acid.xyz"
        reference.write_text("2\nref\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")
        batch_cmd = pka_mod.pka.commands["batch"]

        shared = _build_pka_shared(
            scheme="proton exchange",
            reference=str(reference),
            reference_charge=0,
            reference_multiplicity=1,
        )
        ctx = click.Context(batch_cmd)
        ctx.obj = {
            "pka_shared": shared,
            "filename": str(table),
            "jobrunner": None,
            "project_settings": project_settings,
            "job_settings": None,
            "keywords": {},
        }
        with ctx:
            with pytest.raises(click.UsageError) as exc_info:
                batch_cmd.callback(
                    skip_completed=False, proton_index=None, color_code=None
                )
        assert "-rpi/--reference-proton-index" in str(exc_info.value)

    def test_batch_proton_exchange_cdxml_reference_resolve_failure_raises_usage_error(
        self, tmp_path, monkeypatch
    ):
        """Documents BUGS_FOUND.md #42: resolve_reference_proton raises
        ValueError on failure, but batch() only catches click.UsageError,
        so a real CDXML reference-proton auto-detect failure currently
        propagates uncaught instead of becoming a friendly UsageError."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        table = _build_pka_batch_table(tmp_path)

        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")
        batch_cmd = pka_mod.pka.commands["batch"]

        shared = _build_pka_shared(
            scheme="proton exchange",
            reference="reference.cdxml",
            reference_charge=0,
            reference_multiplicity=1,
        )
        ctx = click.Context(batch_cmd)
        ctx.obj = {
            "pka_shared": shared,
            "filename": str(table),
            "jobrunner": None,
            "project_settings": project_settings,
            "job_settings": None,
            "keywords": {},
        }
        with (
            ctx,
            monkeypatch.context() as m,
        ):
            m.setattr(
                pka_mod.PKaCDXFile,
                "resolve_reference_proton",
                staticmethod(
                    lambda *a, **k: (_ for _ in ()).throw(
                        ValueError("no coloured proton found")
                    )
                ),
            )
            with pytest.raises(ValueError, match="no coloured proton found"):
                batch_cmd.callback(
                    skip_completed=False, proton_index=None, color_code=None
                )

    def test_batch_proton_exchange_full_reference_options_succeeds(
        self, tmp_path, monkeypatch
    ):
        """All reference options supplied: the missing-options check
        passes (covers the `if missing:` False arm) and the first table
        row keeps scheme='proton exchange' while later rows are forced
        to 'direct' (covers both arms of the per-row scheme rewrite)."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        table = _build_pka_batch_table(tmp_path)
        reference = tmp_path / "reference_acid.xyz"
        reference.write_text("2\nref\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")
        batch_cmd = pka_mod.pka.commands["batch"]

        shared = _build_pka_shared(
            scheme="proton exchange",
            reference=str(reference),
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
        )
        ctx = click.Context(batch_cmd)
        ctx.obj = {
            "pka_shared": shared,
            "filename": str(table),
            "jobrunner": None,
            "project_settings": project_settings,
            "job_settings": None,
            "keywords": {},
        }
        with ctx:
            jobs = batch_cmd.callback(
                skip_completed=False, proton_index=None, color_code=None
            )

        assert len(jobs) == 2
        assert jobs[0]._batch_entry["scheme"] == "proton exchange"
        assert jobs[1]._batch_entry["scheme"] == "direct"

    def test_batch_job_settings_none_skips_merge(self, tmp_path, monkeypatch):
        """ctx.obj['job_settings'] falsy should skip the opt_settings.merge
        call rather than erroring, covering batch()'s defensive
        `if job_settings:` False arm."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        table = _build_pka_batch_table(tmp_path)

        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")
        batch_cmd = pka_mod.pka.commands["batch"]

        shared = _build_pka_shared(scheme="direct")
        ctx = click.Context(batch_cmd)
        ctx.obj = {
            "pka_shared": shared,
            "filename": str(table),
            "jobrunner": None,
            "project_settings": project_settings,
            "job_settings": None,
            "keywords": {},
        }
        with ctx:
            jobs = batch_cmd.callback(
                skip_completed=False, proton_index=None, color_code=None
            )

        assert len(jobs) == 2
        assert {job.label for job in jobs} == {"acid1", "acid2"}

    def test_create_pka_jobs_from_molecules_job_settings_none_skips_merge(
        self, tmp_path, monkeypatch
    ):
        """job_settings falsy should skip opt_settings.merge in
        _create_pka_jobs_from_molecules too (its own copy of the same
        defensive guard as batch())."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib
        from types import SimpleNamespace

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")

        water = Molecule(
            symbols=["O", "H"],
            positions=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        water.proton_index = 2

        ctx = SimpleNamespace(
            obj={
                "project_settings": project_settings,
                "job_settings": None,
                "keywords": {},
                "jobrunner": None,
                "filename": "frag.cdxml",
            }
        )
        shared = _build_pka_shared(scheme="direct")

        jobs = pka_mod._create_pka_jobs_from_molecules(
            ctx, [water], shared, False
        )

        assert len(jobs) == 1
        assert jobs[0].label == "frag_frag1_pka"
        assert jobs[0].settings.charge == 0
        assert jobs[0].settings.multiplicity == 1

    def test_submit_empty_molecules_skips_charge_multiplicity_inference(
        self, tmp_path, monkeypatch
    ):
        """ctx.obj['molecules'] falsy (empty list) should skip
        apply_pka_molecule_charge_multiplicity, covering submit()'s
        `if molecules:` False arm. Downstream `molecules[-1]` access
        still requires a non-empty list, so this is expected to fail
        past that point with an IndexError -- the branch under test
        has already run by then."""
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.gaussian import GaussianProjectSettings

        project_settings = GaussianProjectSettings.from_project("test")

        import importlib

        pka_mod = importlib.import_module("chemsmart.cli.gaussian.pka")
        submit_cmd = pka_mod.pka.commands["submit"]

        job_settings = project_settings.opt_settings()
        job_settings.charge = 0
        job_settings.multiplicity = 1

        shared = _build_pka_shared(scheme="direct")
        ctx = click.Context(submit_cmd)
        ctx.obj = {
            "pka_shared": shared,
            "filename": "acid.xyz",
            "jobrunner": None,
            "project_settings": project_settings,
            "job_settings": job_settings,
            "keywords": {"charge", "multiplicity"},
            "molecules": [],
            "molecule_indices": None,
            "label": "acid",
        }
        with ctx:
            with pytest.raises(IndexError):
                submit_cmd.callback(
                    skip_completed=False, proton_index=2, color_code=None
                )


def _invoke_sub_process_pipeline_directly(
    cli_tokens,
    batch_entry=None,
    subcommand_list=None,
    _omit_subcommand_key=False,
    **result_kwargs,
):
    """Directly invoke ``sub``'s result callback with a stubbed
    ``CtxObjArguments.reconstruct_command_line`` output.

    This bypasses the need to build a fully realistic
    ``ctx.obj["subcommand"]`` chain (which only ``MyGroup``/``MyCommand``
    invocation normally populates) so that ``_replace_batch_table_tokens``
    can be driven with exact, hand-picked token lists to reach its less
    common rewrite branches.
    """
    from types import SimpleNamespace
    from unittest.mock import MagicMock, patch

    import click

    import chemsmart.cli.sub as sub_module
    from chemsmart.settings.server import Server

    class _StubArgs:
        def __init__(self, commands, entry_point):
            pass

        def reconstruct_command_line(self):
            return ["placeholder"] + list(cli_tokens)

    fake_server = Server(name="dummy")
    captured = {"cli_args": None}
    fake_server.submit = (
        lambda job, test=False, cli_args=None, **kw: captured.update(
            cli_args=cli_args
        )
    )

    job = SimpleNamespace()
    if batch_entry is not None:
        job._batch_entry = batch_entry

    if subcommand_list is None:
        subcommand_list = [{"name": "sub", "kwargs": {}}]

    kwargs = {"test": False, "print_command": False, "server": "dummy"}
    kwargs.update(result_kwargs)

    with (
        patch.object(sub_module, "CtxObjArguments", _StubArgs),
        patch(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        ),
    ):
        ctx = click.Context(sub_module.sub)
        ctx.obj = {"jobrunner": MagicMock()}
        if not _omit_subcommand_key:
            ctx.obj["subcommand"] = subcommand_list
        with ctx:
            sub_module.sub._result_callback(job, **kwargs)

    return captured["cli_args"]


class TestSubProcessPipelineDirectInvocation:
    """Direct-invocation tests for ``cli/sub.py``'s ``process_pipeline``
    and its nested ``_replace_batch_table_tokens`` helper, covering
    rewrite branches that real pKa-batch CLI runs don't naturally hit."""

    def test_missing_subcommand_key_skips_cleanup(self):
        """``_clean_command`` must no-op (not raise) when ``ctx.obj``
        has no "subcommand" key at all -- ``MyGroup``/``MyCommand``
        always populate it in real invocations, so downstream code
        (``_reconstruct_cli_args``) still assumes it is present and
        raises once cleanup returns; that later ``KeyError`` is exactly
        what proves ``_clean_command`` itself didn't raise first."""
        with pytest.raises(KeyError, match="subcommand"):
            _invoke_sub_process_pipeline_directly(
                ["gaussian", "-p", "test", "opt"],
                subcommand_list=None,
                _omit_subcommand_key=True,
            )

    def test_no_sub_entry_in_subcommand_list_skips_cleanup(self):
        """``_clean_command_list`` finds no entry named "sub" (e.g. if
        invoked from something other than the ``sub`` entry point), so
        it must no-op instead of raising."""
        cli_args = _invoke_sub_process_pipeline_directly(
            ["gaussian", "-p", "test", "opt"],
            subcommand_list=[{"name": "other", "kwargs": {"verbose": {}}}],
        )
        assert cli_args == ["gaussian", "-p", "test", "opt"]

    def test_batch_entry_without_filename_token_leaves_args_unmatched(self):
        """When reconstructed args contain no "-f"/"--filename" token,
        the table-filename-replacement loop must run to completion
        without ever matching (no IndexError, no replacement)."""
        cli_args = _invoke_sub_process_pipeline_directly(
            ["gaussian", "-p", "test", "opt"],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": None,
                "label": None,
            },
        )
        # No "-f" was present to rewrite, and "batch"/"pka" markers are
        # absent too, so the args pass through with only the always-run
        # charge/multiplicity/proton-index rewrites applied (each
        # appended since their markers are absent).
        assert "acid1.xyz" not in cli_args

    def test_batch_scheme_and_label_none_skip_their_rewrite_blocks(self):
        """``scheme``/``label`` are always populated by the real pKa
        job-creation call sites, but the helper still defensively
        guards against them being ``None`` -- exercise that guard
        directly. See BUGS_FOUND.md #54."""
        cli_args = _invoke_sub_process_pipeline_directly(
            ["gaussian", "-p", "test", "-f", "table.csv", "pka", "submit"],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": None,
                "label": None,
            },
        )
        assert "--label" not in cli_args
        assert "--scheme" not in cli_args
        assert "acid1.xyz" in cli_args

    def test_batch_scheme_direct_drops_trailing_reference_flag_without_value(
        self,
    ):
        """``_drop_option`` must handle a reference flag that is the
        very last token (no following value) without raising."""
        cli_args = _invoke_sub_process_pipeline_directly(
            [
                "gaussian",
                "-p",
                "test",
                "-f",
                "table.csv",
                "pka",
                "submit",
                "--reference",
            ],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": "direct",
                "label": "acid1_pka",
            },
        )
        assert "--reference" not in cli_args

    def test_batch_entry_replaces_preexisting_proton_index_tokens(self):
        """If the reconstructed args already contain a stale
        "--proton-index" pair (e.g. from the table's own submit-level
        option), it must be dropped and replaced by the row-specific
        value, not duplicated."""
        cli_args = _invoke_sub_process_pipeline_directly(
            [
                "gaussian",
                "-p",
                "test",
                "-f",
                "table.csv",
                "pka",
                "submit",
                "--proton-index",
                "99",
            ],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": "direct",
                "label": "acid1_pka",
            },
        )
        assert cli_args.count("--proton-index") == 1
        idx = cli_args.index("--proton-index")
        assert cli_args[idx + 1] == "2"

    def test_batch_entry_charge_as_trailing_token_without_value(self):
        """``_set_option`` must handle its long option being the very
        last token in the reconstructed args (no following value).

        "submit" must precede "--charge" here so that the always-run
        proton-index insertion (which targets the position right after
        "submit") lands before "--charge" instead of after it, keeping
        "--charge" the last token by the time ``_set_option`` runs."""
        cli_args = _invoke_sub_process_pipeline_directly(
            [
                "gaussian",
                "-p",
                "test",
                "-f",
                "table.csv",
                "pka",
                "submit",
                "--charge",
            ],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": None,
                "label": None,
            },
        )
        assert cli_args.count("--charge") == 1
        assert cli_args[-1] == "--charge"

    def test_batch_entry_short_form_charge_flag_gets_updated_value(self):
        """``_set_option`` must also resolve the short-form alias
        ("-c") when the long form isn't present."""
        cli_args = _invoke_sub_process_pipeline_directly(
            [
                "gaussian",
                "-p",
                "test",
                "-f",
                "table.csv",
                "-c",
                "99",
                "pka",
            ],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": None,
                "label": None,
            },
        )
        idx = cli_args.index("-c")
        assert cli_args[idx + 1] == "0"

    def test_batch_entry_short_form_charge_as_trailing_token_without_value(
        self,
    ):
        """``_set_option`` must handle the short-form alias being the
        very last token too (no following value), mirroring the
        long-option case above.

        A pre-existing "-m 1" pair is included so the always-run
        multiplicity rewrite updates it in place instead of appending
        a new pair after "-c", which would otherwise disturb the
        "-c is last" setup this test relies on."""
        cli_args = _invoke_sub_process_pipeline_directly(
            [
                "gaussian",
                "-p",
                "test",
                "-f",
                "table.csv",
                "-m",
                "1",
                "submit",
                "-c",
            ],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": None,
                "label": None,
            },
        )
        assert cli_args.count("-c") == 1
        assert cli_args[-1] == "-c"

    def test_batch_entry_without_pka_marker_appends_options_at_end(self):
        """When the "pka" insertion marker is absent from the
        reconstructed args, charge/multiplicity options must be
        appended at the end instead of raising."""
        cli_args = _invoke_sub_process_pipeline_directly(
            ["gaussian", "-p", "test", "-f", "table.csv"],
            batch_entry={
                "filepath": "acid1.xyz",
                "proton_index": 2,
                "charge": 0,
                "multiplicity": 1,
                "scheme": None,
                "label": None,
            },
        )
        assert cli_args[-4:] == [
            "--charge",
            "0",
            "--multiplicity",
            "1",
        ]

    def test_print_command_flag_prints_reconstructed_args(self, capsys):
        _invoke_sub_process_pipeline_directly(
            ["gaussian", "-p", "test", "opt"],
            print_command=True,
        )
        captured = capsys.readouterr()
        assert "gaussian" in captured.out

    def test_non_test_submission_skips_test_warning_log(self, caplog):
        _invoke_sub_process_pipeline_directly(
            ["gaussian", "-p", "test", "opt"],
            test=False,
        )
        assert "Not submitting" not in caplog.text
