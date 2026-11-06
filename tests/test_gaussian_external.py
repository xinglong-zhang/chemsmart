"""
Tests for running Gaussian optimisers on an External level of theory.

Covers:
  chemsmart.io.gaussian.external    file I/O, engines, Hessians, driver
  GaussianExternalCalculator        validation, route value, launcher
  GaussianJobSettings.external      route strings for opt/TS/IRC/sp
  GaussianJobRunner                 writes the launcher next to the .com
  ``gaussian ... external`` CLI     reuses opt/ts/... with the backend
  Real g16 runs                     skipped unless Gaussian is installed
"""

import os
import subprocess
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from ase.units import Bohr, Hartree

from chemsmart.io.gaussian.external import (
    ASEEngine,
    XTBEngine,
    build_engine,
    compute_hessian,
    finite_difference_hessian,
    main,
    parse_xtb_gradient_file,
    parse_xtb_hessian_file,
    read_external_input,
    run,
    write_external_output,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.external import GaussianExternalCalculator
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
)
from tests.helpers import (
    QuadraticEngine,
    run_gaussian_external_cli,
    run_h2_g16_external,
)

# ── File formats ──────────────────────────────────────────────────────────


class TestExternalFileFormats:
    def test_read_input_converts_bohr_to_angstrom(self, au2_external_input):
        request = read_external_input(
            au2_external_input(deriv=2, charge=-1, multiplicity=2)
        )
        assert request.natoms == 2
        assert request.deriv == 2
        assert request.charge == -1
        assert request.multiplicity == 2
        assert list(request.numbers) == [79, 79]
        np.testing.assert_allclose(request.positions[1], [0, 0, 2.88])

    def test_energy_only_output(self, tmp_path):
        out = tmp_path / "out"
        write_external_output(out, -1.5)
        lines = out.read_text().splitlines()
        assert len(lines) == 1
        assert float(lines[0].split()[0]) == pytest.approx(-1.5)

    def test_second_derivative_block_layout(self, tmp_path, harmonic_hessian):
        out = tmp_path / "out"
        gradient = np.arange(6.0).reshape(2, 3)
        write_external_output(out, -1.0, gradient, harmonic_hessian)
        lines = out.read_text().splitlines()

        n = 6
        n_lower = n * (n + 1) // 2
        # energy, 2 gradient rows, 2 polar rows, 3N dipole-derivative
        # rows, then the lower-triangle force constants 3 per line.
        assert len(lines) == 1 + 2 + 2 + n + -(-n_lower // 3)
        fc_lines = lines[1 + 2 + 2 + n :]
        values = [float(v) for line in fc_lines for v in line.split()]
        np.testing.assert_allclose(
            values, harmonic_hessian[np.tril_indices(n)]
        )
        # off-diagonal coupling lands at lower-triangle position (3, 0)
        assert values[6] == pytest.approx(-0.1)


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

    def test_parse_gradient_file(self, xtb_output_directory):
        energy, gradient = parse_xtb_gradient_file(xtb_output_directory, 2)
        assert energy == pytest.approx(-0.9822404321)
        np.testing.assert_allclose(
            gradient, [[0, 0, -0.012345678901234], [0, 0, 0.012345678901234]]
        )

    def test_parse_hessian_file(self, xtb_output_directory):
        hessian = parse_xtb_hessian_file(xtb_output_directory, 2)
        assert hessian.shape == (6, 6)
        assert hessian[1, 0] == pytest.approx(6.0)

    def test_parse_hessian_rejects_wrong_size(self, xtb_output_directory):
        with pytest.raises(ValueError, match="expected 9"):
            parse_xtb_hessian_file(xtb_output_directory, 1)


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

    def test_analytic_mode_requires_analytic_hessian(
        self, harmonic_hessian, harmonic_external_request
    ):
        with pytest.raises(RuntimeError, match="analytic"):
            compute_hessian(
                QuadraticEngine(harmonic_hessian),
                harmonic_external_request,
                "analytic",
            )


# ── Backends and driver ───────────────────────────────────────────────────


class TestBuildEngineAndDriver:
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

    def test_run_writes_hessian_block_for_deriv_2(
        self, tmp_path, au2_external_input
    ):
        out = tmp_path / "out"
        message = run(
            {"backend": "ase", "calculator": "ase.calculators.emt.EMT"},
            au2_external_input(deriv=2),
            out,
        )
        assert "Hessian from fd" in message
        assert len(out.read_text().splitlines()) == 1 + 2 + 2 + 6 + 7

    def test_main_reports_failure_in_message_file(self, tmp_path):
        msg = tmp_path / "msg"
        code = main(
            {"backend": "ase", "calculator": "ase.calculators.emt.EMT"},
            ["R", str(tmp_path / "missing"), "out", str(msg), "f", "m"],
        )
        assert code == 1
        assert "failed" in msg.read_text()

    def test_main_requires_six_gaussian_arguments(self):
        assert main({"backend": "xtb"}, ["only", "three", "args"]) == 1


# ── Backend configuration ─────────────────────────────────────────────────


class TestGaussianExternalCalculator:
    @pytest.mark.parametrize(
        "kwargs, message",
        [
            ({"backend": "orca"}, "Unknown External backend"),
            ({"backend": "ase"}, "dotted path"),
            ({"backend": "xtb", "hessian": "exact"}, "hessian"),
            ({"backend": "xtb", "fd_step": 0}, "fd_step"),
        ],
    )
    def test_validation(self, kwargs, message):
        with pytest.raises(ValueError, match=message):
            GaussianExternalCalculator(**kwargs)

    def test_script_name_depends_only_on_backend_config(
        self, xtb_external_calculator
    ):
        same = GaussianExternalCalculator.from_dict(
            {**xtb_external_calculator.to_dict(), "python_executable": None}
        )
        other = GaussianExternalCalculator(
            backend="xtb", options={"method": "gfn1"}
        )
        name = xtb_external_calculator.script_filename
        assert name.startswith("chemsmart_external_xtb_")
        assert same.script_filename == name
        assert other.script_filename != name

    def test_route_value_quotes_interpreter_and_script(
        self, xtb_external_calculator
    ):
        assert xtb_external_calculator.route_value == (
            '"/opt/conda/bin/python '
            f'{xtb_external_calculator.script_filename}"'
        )

    def test_write_script_is_executable_launcher(
        self, tmp_path, emt_external_calculator
    ):
        path = emt_external_calculator.write_script(tmp_path)
        assert os.access(path, os.X_OK)
        content = open(path).read()
        assert content.startswith(f"#!{sys.executable}\n")
        assert "from chemsmart.io.gaussian.external import main" in content
        assert "'ase.calculators.emt.EMT'" in content

    def test_xtb_executable_comes_from_server_config(self, tmp_path):
        calc = GaussianExternalCalculator(backend="xtb")
        executable = MagicMock()
        executable.get_executable.return_value = "/srv/xtb/bin/xtb"
        server = MagicMock()
        server.name = "cluster"
        with patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            return_value=executable,
        ) as from_servername:
            path = calc.write_script(tmp_path, server=server)
        from_servername.assert_called_once_with("cluster")
        assert "'executable': '/srv/xtb/bin/xtb'" in open(path).read()
        # the launcher name is unchanged by server-resolved defaults
        assert os.path.basename(path) == calc.script_filename

    def test_launcher_answers_gaussian_like_gaussian_does(
        self, tmp_path, emt_external_calculator, au2_external_input
    ):
        """Run the launcher exactly as Gaussian would, requesting a Hessian."""
        script = emt_external_calculator.write_script(tmp_path)
        out, msg = tmp_path / "ext.out", tmp_path / "ext.msg"
        result = subprocess.run(
            [
                sys.executable,
                script,
                "R",
                au2_external_input(deriv=2),
                str(out),
                str(msg),
                str(tmp_path / "ext.fchk"),
                str(tmp_path / "ext.matel"),
            ],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert result.returncode == 0, result.stderr + msg.read_text()
        lines = out.read_text().splitlines()
        assert len(lines) == 1 + 2 + 2 + 6 + 7
        assert abs(float(lines[0].split()[0])) < 10.0  # Hartree, not eV


# ── Settings, writer and runner ───────────────────────────────────────────


class TestExternalSettings:
    @pytest.mark.parametrize(
        "jobtype, expected_prefix",
        [
            ("opt", "# opt freq External="),
            ("ts", "# opt=(ts,calcfc,noeigentest) freq External="),
            ("sp", "# External="),
        ],
    )
    def test_external_replaces_level_of_theory(
        self,
        external_dft_settings,
        xtb_external_calculator,
        jobtype,
        expected_prefix,
    ):
        external_dft_settings.jobtype = jobtype
        settings = external_dft_settings.merge(
            {"external": xtb_external_calculator}, keywords=("external",)
        )
        route = settings.route_string
        assert route == expected_prefix + xtb_external_calculator.route_value
        assert "b3lyp" not in route

    def test_irc_route_uses_external(self, xtb_external_calculator):
        settings = GaussianIRCJobSettings(
            functional="b3lyp",
            basis="def2svp",
            jobtype="ircf",
            external=xtb_external_calculator.to_dict(),
        )
        assert settings.external == xtb_external_calculator
        assert "External=" in settings.route_string
        assert "irc(calcfc" in settings.route_string
        assert "b3lyp" not in settings.route_string

    def test_gaussian_solvation_and_genecp_dropped(
        self, xtb_external_calculator, caplog
    ):
        settings = GaussianJobSettings(
            functional="b3lyp",
            basis="genecp",
            heavy_elements=["Pd"],
            heavy_elements_basis="def2-tzvp",
            light_elements_basis="def2svp",
            solvent_model="smd",
            solvent_id="water",
            jobtype="opt",
            external=xtb_external_calculator,
        )
        assert "scrf" not in settings.route_string
        assert settings.genecp is False
        assert "solvation is ignored" in caplog.text

    def test_copy_preserves_external(self, xtb_external_calculator):
        settings = GaussianJobSettings(external=xtb_external_calculator)
        assert settings.copy().external == xtb_external_calculator

    def test_runner_writes_launcher_next_to_input(
        self,
        tmp_path,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        emt_external_calculator,
    ):
        settings = GaussianJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            external=emt_external_calculator,
        )
        job = GaussianOptJob(
            molecule=Molecule.from_filepath(single_molecule_xyz_file),
            settings=settings,
            label="mol_opt_ase",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)
        gaussian_jobrunner_no_scratch._prerun(job)
        gaussian_jobrunner_no_scratch._write_input(job)

        com = open(gaussian_jobrunner_no_scratch.job_inputfile).read()
        assert f"External={emt_external_calculator.route_value}" in com
        assert (tmp_path / emt_external_calculator.script_filename).is_file()


# ── CLI ───────────────────────────────────────────────────────────────────


class TestExternalCLI:
    def test_ts_with_xtb(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        job_cls = MagicMock(name="GaussianTSJob")
        result, settings = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.ts.GaussianTSJob",
            "-b",
            "xtb",
            "-o",
            "method=gfn2",
            "-o",
            "alpb=water",
            "ts",
            job_cls=job_cls,
        )
        assert result.exit_code == 0, result.output
        assert settings.external == GaussianExternalCalculator(
            backend="xtb", options={"method": "gfn2", "alpb": "water"}
        )
        assert settings.route_string.startswith(
            "# opt=(ts,calcfc,noeigentest)"
        )
        assert job_cls.call_args.kwargs["label"] == "crest_best_ts_xtb"

    def test_opt_with_ase_calculator_and_typed_options(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        result, settings = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            "-b",
            "ase",
            "--calculator",
            "ase.calculators.lj.LennardJones",
            "-o",
            "sigma=0.5",
            "--hessian",
            "fd",
            "opt",
        )
        assert result.exit_code == 0, result.output
        external = settings.external
        assert external.calculator == "ase.calculators.lj.LennardJones"
        assert external.options == {"sigma": 0.5}
        assert external.hessian == "fd"

    def test_malformed_option_is_rejected(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        result, _ = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            "-b",
            "xtb",
            "-o",
            "gfn2",
            "opt",
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "key=value" in result.output

    def test_ase_backend_requires_calculator(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        result, _ = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            "-b",
            "ase",
            "opt",
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "dotted path" in result.output


# ── Real Gaussian runs ────────────────────────────────────────────────────


class TestGaussianExternalRealRun:
    """
    Integration tests that invoke the real g16 binary.

    Each test is skipped automatically when Gaussian is not installed
    (via the ``g16_executable`` fixture).

    ASE's ``EMT`` cannot evaluate H, so H2 uses ASE's ``LennardJones``
    with sigma=0.5 Angstrom, which puts r=0.74 Angstrom in the attractive
    well. LennardJones has no analytic Hessian, so ``calcfc`` exercises
    the finite-difference second-derivative path.
    """

    def test_h2_lj_external_singlepoint(
        self, tmp_path, g16_executable, lj_external_calculator
    ):
        log_text = run_h2_g16_external(
            g16_executable, tmp_path, lj_external_calculator, "", 0.74, "h2_sp"
        )
        assert "Energy=" in log_text

    def test_h2_lj_external_opt(
        self, tmp_path, g16_executable, lj_external_calculator
    ):
        log_text = run_h2_g16_external(
            g16_executable,
            tmp_path,
            lj_external_calculator,
            "opt",
            0.90,
            "h2_opt",
        )
        assert "Stationary point found" in log_text

    def test_h2_lj_external_opt_calcfc(
        self, tmp_path, g16_executable, lj_external_calculator
    ):
        """``calcfc`` makes Gaussian request a Hessian (deriv=2), the same
        request a TS search makes."""
        log_text = run_h2_g16_external(
            g16_executable,
            tmp_path,
            lj_external_calculator,
            "opt=calcfc",
            0.90,
            "h2_calcfc",
        )
        assert "Stationary point found" in log_text
