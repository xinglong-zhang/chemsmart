"""
Tests for the program-neutral external backends in
:mod:`chemsmart.io.external`: ASE and xtb engines, ML-potential presets,
and analytic/finite-difference Hessians. Shared by the Gaussian
(``External=``) and ORCA (``ExtOpt``) interfaces.
"""

import importlib.metadata
import json
import re
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import click
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from ase.units import Bohr, Hartree

from chemsmart.cli.external import reject_external_level_of_theory
from chemsmart.io.external import (
    ENTRY_POINT_GROUP,
    ASEEngine,
    ExternalEngine,
    ExternalProtocol,
    XTBEngine,
    available_backends,
    build_engine,
    compute_hessian,
    evaluate,
    finite_difference_hessian,
    get_backend,
    read_xtb_gradient,
    read_xtb_hessian,
    register_backend,
    validate_backend,
)
from chemsmart.io.external.runlog import EXTERNAL_LOG_ENV, gradient_statistics
from chemsmart.utils.repattern import external_log_step_pattern
from tests.helpers import (
    JSONExternalProtocol,
    QuadraticEngine,
    read_external_log_records,
    run_json_external_request,
    write_json_external_request,
)

# ── Engines ───────────────────────────────────────────────────────────────


class TestASEEngine:
    def test_energy_and_gradient_units_and_sign(self):
        positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.88]]
        reference = Atoms("Au2", positions=positions)
        reference.calc = EMT()
        engine = ASEEngine(lambda charge, multiplicity: EMT())

        energy, gradient = engine.energy_and_gradient(
            [79, 79], positions, 0, 1
        )

        assert energy == pytest.approx(
            reference.get_potential_energy() / Hartree
        )
        np.testing.assert_allclose(
            gradient, -reference.get_forces() * Bohr / Hartree
        )

    def test_charge_and_spin_reach_calculator_and_atoms(self):
        seen = {}

        def factory(charge, multiplicity):
            seen["args"] = (charge, multiplicity)
            calc = EMT()
            original = calc.calculate

            def calculate(atoms=None, *args, **kwargs):
                seen["info"] = dict(atoms.info)
                return original(atoms, *args, **kwargs)

            calc.calculate = calculate
            return calc

        ASEEngine(factory).energy_and_gradient(
            [79, 79], [[0, 0, 0], [0, 0, 2.88]], -1, 2
        )
        assert seen["args"] == (-1, 2)
        assert seen["info"] == {"charge": -1, "spin": 2}

    def test_no_get_hessian_means_no_analytic_hessian(self):
        engine = ASEEngine(lambda charge, multiplicity: EMT())
        assert (
            engine.analytic_hessian([79, 79], [[0, 0, 0], [0, 0, 2.88]], 0, 1)
            is None
        )

    def test_get_hessian_converted_to_atomic_units(self):
        hessian_ev = np.arange(36.0).reshape(6, 6)

        def factory(charge, multiplicity):
            calc = EMT()
            # flat output must be reshaped to (3N, 3N)
            calc.get_hessian = lambda atoms: hessian_ev.ravel()
            return calc

        hessian = ASEEngine(factory).analytic_hessian(
            [79, 79], [[0, 0, 0], [0, 0, 2.88]], 0, 1
        )
        np.testing.assert_allclose(hessian, hessian_ev * Bohr**2 / Hartree)

    def test_failing_get_hessian_means_no_analytic_hessian(self, caplog):
        def factory(charge, multiplicity):
            calc = EMT()
            calc.get_hessian = MagicMock(side_effect=RuntimeError("no model"))
            return calc

        hessian = ASEEngine(factory).analytic_hessian(
            [79, 79], [[0, 0, 0], [0, 0, 2.88]], 0, 1
        )
        assert hessian is None
        assert "no model" in caplog.text

    def test_base_engine_contract(self):
        engine = ExternalEngine()
        with pytest.raises(NotImplementedError):
            engine.energy_and_gradient([1], [[0, 0, 0]], 0, 1)
        assert engine.analytic_hessian([1], [[0, 0, 0]], 0, 1) is None


class TestXTBEngine:
    def test_command_carries_charge_uhf_method_and_solvent(self):
        engine = XTBEngine(
            executable="/opt/xtb/bin/xtb", method="GFN2-xTB", alpb="water"
        )
        command = engine._command("--grad", charge=-1, multiplicity=2)
        assert command == [
            "/opt/xtb/bin/xtb",
            "coord.xyz",
            "--grad",
            "--chrg",
            "-1",
            "--uhf",
            "1",
            "--gfn",
            "2",
            "--alpb",
            "water",
        ]

    @pytest.mark.parametrize(
        "kwargs",
        [{"method": "pm7"}, {"alpb": "water", "gbsa": "water"}],
    )
    def test_invalid_configuration_rejected(self, kwargs):
        with pytest.raises(ValueError):
            XTBEngine(**kwargs)

    def test_gbsa_and_string_args_in_command(self):
        engine = XTBEngine(gbsa="thf", args="--acc 0.1 --etemp 500")
        command = engine._command("--hess", charge=0, multiplicity=1)
        assert command[-6:] == [
            "--gbsa",
            "thf",
            "--acc",
            "0.1",
            "--etemp",
            "500",
        ]
        assert "--alpb" not in command

    def test_energy_and_gradient_runs_xtb(self, fake_xtb_executable):
        executable, log = fake_xtb_executable
        energy, gradient = XTBEngine(
            executable=executable
        ).energy_and_gradient([1, 1], [[0, 0, 0], [0, 0, 0.74]], 0, 1)
        assert energy == pytest.approx(-0.9822404321)
        assert gradient.shape == (2, 3)
        calls = log.read_text()
        assert "coord.xyz --grad --chrg 0 --uhf 0 --gfn 2" in calls
        assert "H " in calls and "0.7400000000" in calls
        # temporary xtb_* working folders are removed afterwards
        assert not list(log.parent.glob("xtb_*"))

    def test_analytic_hessian_runs_xtb(self, fake_xtb_executable):
        executable, log = fake_xtb_executable
        hessian = XTBEngine(executable=executable).analytic_hessian(
            [1, 1], [[0, 0, 0], [0, 0, 0.74]], 0, 1
        )
        assert hessian.shape == (6, 6)
        assert "--hess" in log.read_text()

    def test_xtb_failure_raises_with_output_tail(self, fake_xtb_executable):
        executable, _ = fake_xtb_executable
        engine = XTBEngine(executable=executable, args=["--fail"])
        with pytest.raises(RuntimeError, match=r"exit 3(.|\n)*xtb: boom"):
            engine.energy_and_gradient([1, 1], [[0, 0, 0], [0, 0, 0.74]], 0, 1)

    def test_xtb_hessian_failure_returns_none(
        self, fake_xtb_executable, caplog
    ):
        executable, _ = fake_xtb_executable
        engine = XTBEngine(executable=executable, args=["--fail"])
        hessian = engine.analytic_hessian(
            [1, 1], [[0, 0, 0], [0, 0, 0.74]], 0, 1
        )
        assert hessian is None
        assert "xtb Hessian unavailable" in caplog.text

    def test_read_gradient_file(self, xtb_output_directory):
        energy, gradient = read_xtb_gradient(xtb_output_directory, 2)
        assert energy == pytest.approx(-0.9822404321)
        np.testing.assert_allclose(
            gradient, [[0, 0, -0.012345678901234], [0, 0, 0.012345678901234]]
        )

    def test_read_hessian_file(self, xtb_output_directory):
        hessian = read_xtb_hessian(xtb_output_directory, 2)
        assert hessian.shape == (6, 6)
        assert hessian[1, 0] == pytest.approx(6.0)

    def test_read_hessian_rejects_wrong_size(self, xtb_output_directory):
        with pytest.raises(ValueError, match=r"expected \(3, 3\)"):
            read_xtb_hessian(xtb_output_directory, 1)

    def test_read_gradient_rejects_truncated_file(self, tmp_path):
        (tmp_path / "gradient").write_text(
            "$grad\n"
            "  cycle = 1    SCF energy =    -1.0   |dE/xyz| =  0.0\n"
            "  0.0  0.0  0.0  h\n"
            "  0.0  0.0  1.4  h\n"
            "  0.0  0.0  -1.0E-02\n"
            "$end\n"
        )
        with pytest.raises(ValueError, match="gradient shape"):
            read_xtb_gradient(str(tmp_path), 2)

    def test_read_gradient_requires_energy_and_gradient(self, tmp_path):
        (tmp_path / "gradient").write_text("$grad\n$end\n")
        with pytest.raises(ValueError, match="No energy/gradient"):
            read_xtb_gradient(str(tmp_path), 2)


# ── Hessians ──────────────────────────────────────────────────────────────


class TestHessian:
    def test_finite_difference_recovers_harmonic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        engine = QuadraticEngine(harmonic_hessian)
        hessian = finite_difference_hessian(engine, harmonic_external_request)
        np.testing.assert_allclose(hessian, harmonic_hessian, atol=1e-8)
        assert engine.gradient_calls == 2 * 6

    def test_auto_uses_consistent_analytic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        engine = QuadraticEngine(
            harmonic_hessian, analytic=harmonic_hessian.copy()
        )
        hessian, source = compute_hessian(engine, harmonic_external_request)
        assert source == "analytic"
        np.testing.assert_allclose(hessian, harmonic_hessian)
        assert engine.gradient_calls == 2  # one spot-check column

    def test_auto_rejects_wrong_sign_analytic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        engine = QuadraticEngine(harmonic_hessian, analytic=-harmonic_hessian)
        hessian, source = compute_hessian(engine, harmonic_external_request)
        assert source == "fd"
        np.testing.assert_allclose(hessian, harmonic_hessian, atol=1e-8)

    def test_auto_falls_back_without_analytic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        _, source = compute_hessian(
            QuadraticEngine(harmonic_hessian), harmonic_external_request
        )
        assert source == "fd"

    def test_fd_mode_ignores_analytic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        engine = QuadraticEngine(harmonic_hessian, analytic=-harmonic_hessian)
        _, source = compute_hessian(
            engine, harmonic_external_request, mode="fd"
        )
        assert source == "fd"

    def test_analytic_mode_trusts_engine_without_spot_check(
        self, harmonic_hessian, harmonic_external_request
    ):
        asymmetric = harmonic_hessian.copy()
        asymmetric[0, 3] = 0.0  # symmetrised to -0.05
        engine = QuadraticEngine(harmonic_hessian, analytic=asymmetric)
        hessian, source = compute_hessian(
            engine, harmonic_external_request, mode="analytic"
        )
        assert source == "analytic"
        assert hessian[0, 3] == hessian[3, 0] == pytest.approx(-0.05)
        assert engine.gradient_calls == 0

    def test_unknown_mode_rejected(
        self, harmonic_hessian, harmonic_external_request
    ):
        with pytest.raises(ValueError, match="hessian mode"):
            compute_hessian(
                QuadraticEngine(harmonic_hessian),
                harmonic_external_request,
                mode="exact",
            )

    def test_analytic_mode_requires_analytic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        with pytest.raises(RuntimeError, match="analytic"):
            compute_hessian(
                QuadraticEngine(harmonic_hessian),
                harmonic_external_request,
                "analytic",
            )


# ── Backends ──────────────────────────────────────────────────────────────


class TestBuildEngine:
    def test_ase_backend_imports_dotted_path(self):
        engine = build_engine(
            {"backend": "ase", "calculator": "ase.calculators.emt.EMT"}
        )
        assert isinstance(engine, ASEEngine)

    def test_xtb_backend_passes_options(self):
        engine = build_engine(
            {"backend": "xtb", "options": {"method": "gfn1", "gbsa": "thf"}}
        )
        assert isinstance(engine, XTBEngine)
        assert (engine.method, engine.gbsa) == ("gfn1", "thf")

    def test_mace_off_preset_defaults_to_float64(self):
        engine = build_engine({"backend": "mace-off", "options": {}})
        fake_mace = MagicMock()
        with patch.dict(sys.modules, {"mace": fake_mace}):
            engine.factory(0, 1)
        fake_mace.calculators.mace_off.assert_called_once_with(
            model="medium", device="cpu", default_dtype="float64"
        )

    def test_aimnet2_preset_passes_charge_and_multiplicity(self):
        engine = build_engine({"backend": "aimnet2", "options": {}})
        fake = MagicMock()
        with patch.dict(
            sys.modules, {"aimnet": fake, "aimnet.calculators": fake}
        ):
            engine.factory(1, 2)
        fake.AIMNet2ASE.assert_called_once_with("aimnet2", charge=1, mult=2)

    def test_uma_preset_builds_fairchem_calculator(self):
        engine = build_engine(
            {"backend": "uma", "options": {"device": "cuda", "seed": 1}}
        )
        fake = MagicMock()
        with patch.dict(
            sys.modules, {"fairchem": fake, "fairchem.core": fake}
        ):
            calc = engine.factory(0, 1)
        fake.pretrained_mlip.get_predict_unit.assert_called_once_with(
            "uma-s-1", device="cuda"
        )
        fake.FAIRChemCalculator.assert_called_once_with(
            fake.pretrained_mlip.get_predict_unit.return_value,
            task_name="omol",
            seed=1,
        )
        assert calc is fake.FAIRChemCalculator.return_value

    @pytest.mark.parametrize(
        "config, message",
        [
            ({"backend": "orca"}, "Unknown External backend"),
            ({"backend": "ase", "calculator": "EMT"}, "dotted path"),
        ],
    )
    def test_bad_backend_config_rejected(self, config, message):
        with pytest.raises(ValueError, match=message):
            build_engine(config)


# ── Backend registry ──────────────────────────────────────────────────────


class TestBackendRegistry:
    def test_builtin_backends_registered(self):
        assert {"ase", "mace-off", "mace-mp", "aimnet2", "uma", "xtb"} <= set(
            available_backends()
        )

    def test_lookup_is_case_insensitive(self):
        backend = get_backend("XTB")
        assert backend.name == "xtb"
        assert backend.description

    def test_unknown_backend_lists_choices(self):
        with pytest.raises(ValueError, match="choose from .*xtb"):
            get_backend("orca")

    def test_registered_backend_builds_engines(self, harmonic_backend):
        engine = build_engine({"backend": harmonic_backend})
        assert isinstance(engine, QuadraticEngine)

    def test_duplicate_name_rejected_unless_replaced(
        self, isolated_backend_registry
    ):
        with pytest.raises(ValueError, match="already registered"):
            register_backend("xtb")(lambda options, calculator: None)

        @register_backend("xtb", replace=True, description="replacement")
        def build(options, calculator):
            return "replaced"

        assert build_engine({"backend": "xtb"}) == "replaced"

    def test_validate_hook_runs_without_building(self):
        with pytest.raises(ValueError, match="Unknown xtb method"):
            validate_backend("xtb", {"method": "pm7"})
        assert validate_backend("xtb", {"method": "gfn1"}).name == "xtb"

    def test_required_calculator_checked(self):
        with pytest.raises(ValueError, match="needs a calculator"):
            validate_backend("ase")

    def test_entry_point_plugins_load_once(
        self, isolated_backend_registry, monkeypatch, caplog
    ):
        loads = []

        def register_plugin():
            loads.append("good")
            register_backend("plugin-mlp")(lambda options, calculator: None)

        def broken():
            raise ImportError("missing dependency")

        entry_points = [
            SimpleNamespace(
                name="good", value="pkg.good", load=register_plugin
            ),
            SimpleNamespace(name="bad", value="pkg.bad", load=broken),
        ]
        groups = []

        def fake_entry_points(group):
            groups.append(group)
            return entry_points

        monkeypatch.setattr(
            importlib.metadata, "entry_points", fake_entry_points
        )
        monkeypatch.setattr(
            isolated_backend_registry, "_plugins_loaded", False
        )

        assert "plugin-mlp" in available_backends()
        assert get_backend("plugin-mlp").name == "plugin-mlp"
        assert loads == ["good"]  # loaded once, on first use
        assert groups == [ENTRY_POINT_GROUP]
        assert "missing dependency" in caplog.text


# ── Host protocol base ────────────────────────────────────────────────────


class TestExternalProtocol:
    def test_new_host_needs_only_file_exchange(
        self, tmp_path, harmonic_backend, harmonic_hessian, capsys
    ):
        request = write_json_external_request(
            tmp_path / "request.json", [1, 1], [[0, 0, 0], [0, 0, 1.0]]
        )
        result_file = tmp_path / "result.json"
        code = JSONExternalProtocol().main(
            {"backend": harmonic_backend}, [request, str(result_file)]
        )
        assert code == 0
        result = json.loads(result_file.read_text())
        assert np.shape(result["gradient"]) == (2, 3)
        assert result["hessian"] is None
        assert capsys.readouterr().out.startswith("harmonic: E = ")

    def test_hessian_requests_get_a_hessian(
        self, tmp_path, harmonic_backend, harmonic_hessian
    ):
        request = write_json_external_request(
            tmp_path / "request.json",
            [1, 1],
            [[0, 0, 0], [0, 0, 1.0]],
            deriv=2,
        )
        message = JSONExternalProtocol().run(
            {"backend": harmonic_backend, "hessian": "analytic"},
            (request, str(tmp_path / "result.json")),
        )
        assert message.endswith("; Hessian from analytic")
        result = json.loads((tmp_path / "result.json").read_text())
        np.testing.assert_allclose(result["hessian"], harmonic_hessian)

    def test_energy_only_request_has_no_gradient(
        self, harmonic_hessian, harmonic_external_request
    ):
        harmonic_external_request.deriv = 0
        result = evaluate(
            QuadraticEngine(harmonic_hessian), harmonic_external_request
        )
        assert result.gradient is None
        assert result.summary("x") == f"x: E = {result.energy:.10f} Eh"

    def test_failures_reported_with_exit_code(self, tmp_path, capsys):
        code = JSONExternalProtocol().main(
            {"backend": "xtb"},
            [str(tmp_path / "missing.json"), str(tmp_path / "out.json")],
        )
        assert code == 1
        err = capsys.readouterr().err
        assert "CHEMSMART external calculation failed" in err
        assert "missing.json" in err

    def test_usage_errors(self, capsys):
        assert JSONExternalProtocol().main({"backend": "xtb"}, ["one"]) == 1
        assert "Usage" in capsys.readouterr().err

    def test_base_methods_are_abstract(self):
        protocol = ExternalProtocol()
        for method, args in (
            (protocol.parse_arguments, ([],)),
            (protocol.read_request, (None,)),
            (protocol.write_result, (None, None, None)),
        ):
            with pytest.raises(NotImplementedError):
                method(*args)


# ── Per-job external log ──────────────────────────────────────────────────


class TestExternalRunLog:
    def test_nothing_logged_without_environment(
        self, tmp_path, harmonic_backend, monkeypatch
    ):
        monkeypatch.delenv(EXTERNAL_LOG_ENV, raising=False)
        assert (
            run_json_external_request(tmp_path, {"backend": harmonic_backend})
            == 0
        )
        assert not list(tmp_path.glob("*.log"))

    def test_one_record_per_request(
        self, tmp_path, external_run_log, harmonic_backend
    ):
        config = {"backend": harmonic_backend, "hessian": "analytic"}
        assert run_json_external_request(tmp_path, config) == 0
        assert run_json_external_request(tmp_path, config, deriv=2) == 0

        text = external_run_log.read_text()
        assert text.startswith(
            "# CHEMSMART external log: host JSON, backend harmonic\n"
        )
        assert text.count("# config:") == 1  # header written once
        first, second = read_external_log_records(external_run_log)
        assert " step    1  deriv 1  E = " in first
        assert "hessian -" in first
        assert " step    2  deriv 2  " in second
        assert "hessian analytic" in second
        assert re.search(r"\d+\.\d\d s$", first)

    def test_gradient_statistics_logged(
        self, tmp_path, external_run_log, harmonic_backend, harmonic_hessian
    ):
        run_json_external_request(tmp_path, {"backend": harmonic_backend})
        result = json.loads((tmp_path / "request_result.json").read_text())
        gradient = np.asarray(result["gradient"]).ravel()
        (record,) = read_external_log_records(external_run_log)
        assert f"max|g| = {np.max(np.abs(gradient)):.3e}" in record
        assert f"rms|g| = {np.sqrt(np.mean(gradient**2)):.3e}" in record
        assert f"E = {result['energy']:.12f} Eh" in record

    def test_backend_warnings_logged_under_their_request(
        self, tmp_path, external_run_log, inconsistent_hessian_backend
    ):
        run_json_external_request(
            tmp_path, {"backend": inconsistent_hessian_backend}, deriv=2
        )
        lines = external_run_log.read_text().splitlines()
        record = next(
            i for i, line in enumerate(lines) if " step    1 " in line
        )
        assert "hessian fd" in lines[record]
        assert lines[record + 1].startswith("    warning: Analytic Hessian")

    def test_failures_logged_and_numbered(
        self, tmp_path, external_run_log, harmonic_backend
    ):
        run_json_external_request(tmp_path, {"backend": harmonic_backend})
        code = JSONExternalProtocol().main(
            {"backend": harmonic_backend},
            [str(tmp_path / "missing.json"), str(tmp_path / "out.json")],
        )
        assert code == 1
        records = read_external_log_records(external_run_log)
        assert " step    2  FAILED after " in records[1]
        assert "FileNotFoundError" in records[1]
        assert "    Traceback (most recent call last):" in (
            external_run_log.read_text()
        )

    def test_backend_files_kept_per_step_not_in_log(
        self, tmp_path, external_run_log, fake_xtb_executable
    ):
        executable, _ = fake_xtb_executable
        config = {
            "backend": "xtb",
            "options": {"executable": executable},
            "keep_backend_output": True,
        }
        assert run_json_external_request(tmp_path, config) == 0
        assert run_json_external_request(tmp_path, config, name="second") == 0

        records = read_external_log_records(external_run_log)
        folders = [re.search(r"output (\S+)$", r)[1] for r in records]
        assert re.fullmatch(r"job_external/step_0001_\d{8}-\d{6}", folders[0])
        assert folders[1].startswith("job_external/step_0002_")
        grad = tmp_path / folders[0] / "01_grad"
        assert "fake xtb: normal termination" in (grad / "xtb.out").read_text()
        assert "--grad" in (grad / "command.txt").read_text()
        assert {"coord.xyz", "gradient"} <= {p.name for p in grad.iterdir()}
        # the log holds records only, never backend output
        assert "fake xtb" not in external_run_log.read_text()
        assert not list(tmp_path.glob("xtb_*"))

    def test_backend_files_removed_unless_requested(
        self, tmp_path, external_run_log, fake_xtb_executable
    ):
        executable, _ = fake_xtb_executable
        config = {"backend": "xtb", "options": {"executable": executable}}
        assert run_json_external_request(tmp_path, config) == 0
        (record,) = read_external_log_records(external_run_log)
        assert " output " not in record
        assert not (tmp_path / "job_external").exists()
        assert not list(tmp_path.glob("xtb_*"))

    def test_failed_request_keeps_its_backend_files(
        self, tmp_path, external_run_log, fake_xtb_executable
    ):
        executable, _ = fake_xtb_executable
        config = {
            "backend": "xtb",
            "options": {"executable": executable, "args": ["--fail"]},
            "keep_backend_output": True,
        }
        assert run_json_external_request(tmp_path, config) == 1
        (record,) = read_external_log_records(external_run_log)
        assert " step    1  FAILED after " in record
        folder = re.search(r"output (\S+)$", record)[1]
        xtb_out = tmp_path / folder / "01_grad" / "xtb.out"
        assert "xtb: boom" in xtb_out.read_text()

    def test_log_problems_never_fail_the_request(
        self, tmp_path, harmonic_backend, monkeypatch, capsys
    ):
        blocker = tmp_path / "not_a_directory"
        blocker.write_text("")
        monkeypatch.setenv(EXTERNAL_LOG_ENV, str(blocker / "job.log"))
        assert (
            run_json_external_request(tmp_path, {"backend": harmonic_backend})
            == 0
        )
        assert "CHEMSMART external log not written" in capsys.readouterr().err
        assert (tmp_path / "request_result.json").is_file()

    def test_kept_files_removed_when_log_unwritable(
        self, tmp_path, fake_xtb_executable, monkeypatch, capsys
    ):
        blocker = tmp_path / "not_a_directory"
        blocker.write_text("")
        monkeypatch.setenv(EXTERNAL_LOG_ENV, str(blocker / "job.log"))
        executable, _ = fake_xtb_executable
        config = {
            "backend": "xtb",
            "options": {"executable": executable},
            "keep_backend_output": True,
        }
        assert run_json_external_request(tmp_path, config) == 0
        assert "log not written" in capsys.readouterr().err
        assert not list(tmp_path.glob("xtb_*"))

    def test_concurrent_requests_get_unique_steps(
        self, tmp_path, external_run_log, harmonic_backend
    ):
        """ORCA can evaluate numerical-Hessian displacements in parallel;
        the file lock keeps records whole and step numbers unique."""
        from concurrent.futures import ThreadPoolExecutor

        config = {"backend": harmonic_backend}
        with ThreadPoolExecutor(max_workers=8) as pool:
            codes = list(
                pool.map(
                    lambda i: run_json_external_request(
                        tmp_path, config, name=f"request{i}"
                    ),
                    range(16),
                )
            )
        assert codes == [0] * 16
        steps = sorted(
            int(re.match(external_log_step_pattern, line)[1])
            for line in read_external_log_records(external_run_log)
        )
        assert steps == list(range(1, 17))

    def test_gradient_statistics_without_gradient(self):
        assert gradient_statistics(None) == (None, None)
        assert gradient_statistics([]) == (None, None)

    def test_logs_without_file_locking(
        self, tmp_path, external_run_log, harmonic_backend, monkeypatch
    ):
        """Platforms without fcntl (Windows) still get records."""
        from chemsmart.io.external import runlog

        monkeypatch.setattr(runlog, "fcntl", None)
        assert (
            run_json_external_request(tmp_path, {"backend": harmonic_backend})
            == 0
        )
        assert len(read_external_log_records(external_run_log)) == 1


# ── CLI guards ────────────────────────────────────────────────────────────


class TestRejectExternalLevelOfTheory:
    def test_external_parent_rejected(self, xtb_orca_external_calculator):
        settings = SimpleNamespace(external=xtb_orca_external_calculator)
        with pytest.raises(click.UsageError, match="backend 'xtb'"):
            reject_external_level_of_theory(settings, "qmmm")

    @pytest.mark.parametrize(
        "parent", [None, SimpleNamespace(external=None), MagicMock()]
    )
    def test_other_parents_accepted(self, parent):
        """No parent, a host level of theory, or a mock (whose ``external``
        attribute is a truthy MagicMock) is not an external backend."""
        reject_external_level_of_theory(parent, "qmmm")
