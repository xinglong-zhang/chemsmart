"""
Tests for running Gaussian optimisers on an External level of theory.

Covers:
  chemsmart.io.gaussian.external    Gaussian External file I/O and driver
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

from chemsmart.io.external.runlog import EXTERNAL_LOG_ENV
from chemsmart.io.gaussian.external import (
    main,
    read_external_input,
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
    run_gaussian_external_cli,
    run_h2_g16_external,
    write_external_job_outputs,
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

    def test_gradient_without_hessian_has_no_second_derivative_block(
        self, tmp_path
    ):
        out = tmp_path / "out"
        write_external_output(out, -1.0, np.ones((2, 3)), dipole=[1, 2, 3])
        lines = out.read_text().splitlines()
        assert len(lines) == 1 + 2
        assert [float(v) for v in lines[0].split()] == [-1.0, 1.0, 2.0, 3.0]


# ── Driver ────────────────────────────────────────────────────────────────


class TestGaussianExternalDriver:
    def test_energy_only_for_deriv_0(
        self, emt_external_calculator, tmp_path, au2_external_input
    ):
        out, msg = tmp_path / "out", tmp_path / "msg"
        code = main(
            emt_external_calculator.config,
            ["R", au2_external_input(deriv=0), str(out), str(msg), "f", "m"],
        )
        assert code == 0
        assert "Hessian" not in msg.read_text()
        assert len(out.read_text().splitlines()) == 1

    def test_main_writes_status_to_message_file(
        self, emt_external_calculator, tmp_path, au2_external_input
    ):
        out, msg = tmp_path / "out", tmp_path / "msg"
        code = main(
            emt_external_calculator.config,
            [
                "--extra",
                "R",
                au2_external_input(),
                str(out),
                str(msg),
                "f",
                "m",
            ],
        )
        assert code == 0
        assert msg.read_text().startswith("ase: E = ")
        assert len(out.read_text().splitlines()) == 1 + 2

    def test_hessian_block_for_deriv_2(
        self, emt_external_calculator, tmp_path, au2_external_input
    ):
        out, msg = tmp_path / "out", tmp_path / "msg"
        code = main(
            emt_external_calculator.config,
            ["R", au2_external_input(deriv=2), str(out), str(msg), "f", "m"],
        )
        assert code == 0
        assert "Hessian from fd" in msg.read_text()
        assert len(out.read_text().splitlines()) == 1 + 2 + 2 + 6 + 7

    def test_main_reports_failure_in_message_file(
        self, emt_external_calculator, tmp_path
    ):
        msg = tmp_path / "msg"
        code = main(
            emt_external_calculator.config,
            ["R", str(tmp_path / "missing"), "out", str(msg), "f", "m"],
        )
        assert code == 1
        assert "failed" in msg.read_text()

    def test_main_requires_six_gaussian_arguments(self, capsys):
        assert main({"backend": "xtb"}, ["only", "three", "args"]) == 1
        assert "Usage" in capsys.readouterr().err


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

    def test_missing_server_xtb_config_leaves_xtb_on_path(self, tmp_path):
        calc = GaussianExternalCalculator(backend="xtb")
        server = MagicMock()
        server.name = "laptop"
        with patch(
            "chemsmart.settings.executable.XTBExecutable.from_servername",
            side_effect=FileNotFoundError("no xtb section"),
        ):
            path = calc.write_script(tmp_path, server=server)
        assert "'executable'" not in open(path).read()

    def test_not_equal_to_other_types(self, xtb_external_calculator):
        assert xtb_external_calculator != xtb_external_calculator.to_dict()
        assert "GaussianExternalCalculator(" in repr(xtb_external_calculator)

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

    def test_forces_and_additional_route_parameters_kept(
        self, xtb_external_calculator
    ):
        settings = GaussianJobSettings(
            jobtype="sp",
            forces=True,
            additional_route_parameters="nosymm",
            external=xtb_external_calculator,
        )
        assert settings.route_string == (
            f"# External={xtb_external_calculator.route_value} force nosymm"
        )

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

    def test_rerun_starts_external_log_afresh(
        self,
        tmp_path,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        emt_external_calculator,
    ):
        """Without scratch the log lives in the job folder, so a rerun
        would otherwise append to the previous run's log."""
        runner = gaussian_jobrunner_no_scratch
        job = GaussianOptJob(
            molecule=Molecule.from_filepath(single_molecule_xyz_file),
            settings=GaussianJobSettings(
                charge=0, multiplicity=1, external=emt_external_calculator
            ),
            label="mol_opt_ase",
            jobrunner=runner,
        )
        job.folder = str(tmp_path)
        runner._prerun(job)
        log, outputs = write_external_job_outputs(tmp_path, job.label)

        runner._set_aside_previous_external_outputs(job)

        assert not log.exists() and not outputs.exists()
        (backup,) = job._previous_backup_folders()
        assert os.path.isfile(os.path.join(backup, log.name))
        assert os.path.isdir(os.path.join(backup, outputs.name))

    def test_host_level_jobs_leave_same_named_files_alone(
        self, tmp_path, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        runner = gaussian_jobrunner_no_scratch
        job = GaussianOptJob(
            molecule=Molecule.from_filepath(single_molecule_xyz_file),
            settings=GaussianJobSettings(
                functional="b3lyp", basis="def2svp", charge=0, multiplicity=1
            ),
            label="mol_opt",
            jobrunner=runner,
        )
        job.folder = str(tmp_path)
        runner._prerun(job)
        log, _ = write_external_job_outputs(tmp_path, job.label)

        runner._set_aside_previous_external_outputs(job)

        assert log.exists()
        assert job._previous_backup_folders() == []

    def test_runner_points_launcher_at_job_log(
        self,
        tmp_path,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        emt_external_calculator,
        monkeypatch,
    ):
        runner = gaussian_jobrunner_no_scratch
        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        external_job = GaussianOptJob(
            molecule=molecule,
            settings=GaussianJobSettings(
                charge=0, multiplicity=1, external=emt_external_calculator
            ),
            label="mol_opt_ase",
            jobrunner=runner,
        )
        external_job.folder = str(tmp_path)
        runner._prerun(external_job)
        env = runner._update_os_environ(external_job)
        assert env[EXTERNAL_LOG_ENV] == str(
            tmp_path / f"{external_job.label}_external.log"
        )

        # a DFT job never inherits a log path from the environment
        monkeypatch.setenv(EXTERNAL_LOG_ENV, "/elsewhere/other_external.log")
        dft_job = GaussianOptJob(
            molecule=molecule,
            settings=GaussianJobSettings(
                functional="b3lyp", basis="def2svp", charge=0, multiplicity=1
            ),
            label="mol_opt",
            jobrunner=runner,
        )
        dft_job.folder = str(tmp_path)
        runner._prerun(dft_job)
        assert EXTERNAL_LOG_ENV not in runner._update_os_environ(dft_job)


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

    def test_user_label_is_kept(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        job_cls = MagicMock(name="GaussianOptJob")
        result, _ = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            "-b",
            "xtb",
            "opt",
            gaussian_args=("-l", "my_label"),
            job_cls=job_cls,
        )
        assert result.exit_code == 0, result.output
        assert job_cls.call_args.kwargs["label"] == "my_label"

    def test_sp_label_has_no_solvent_tag(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        job_cls = MagicMock(name="GaussianSinglePointJob")
        result, _ = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
            "-b",
            "xtb",
            "sp",
            job_cls=job_cls,
        )
        assert result.exit_code == 0, result.output
        assert job_cls.call_args.kwargs["label"] == "crest_best_sp_xtb"

    def test_registered_backend_accepted_without_cli_changes(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        harmonic_backend,
    ):
        result, settings = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            "-b",
            harmonic_backend,
            "opt",
        )
        assert result.exit_code == 0, result.output
        assert settings.external.backend == "harmonic"

    @pytest.mark.parametrize("parent", ["opt", "ts"])
    def test_qmmm_under_external_is_rejected(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch, parent
    ):
        """``qmmm`` is reachable through the wrapped opt/ts groups but would
        drop the backend (an ONIOM DFT/MM job labelled *_xtb_qmmm)."""
        job_cls = MagicMock(name="GaussianQMMMJob")
        result, settings = run_gaussian_external_cli(
            single_molecule_xyz_file,
            gaussian_jobrunner_no_scratch,
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob",
            "-b",
            "xtb",
            parent,
            "qmmm",
            "-hx",
            "b3lyp",
            "-hb",
            "6-31g",
            "-lff",
            "uff",
            "-ha",
            "1-5",
            job_cls=job_cls,
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "'qmmm' cannot be combined with an external" in result.output
        assert settings is None
        job_cls.assert_not_called()

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
