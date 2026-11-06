"""
Tests for running ORCA optimisers on an external level of theory
(``! ExtOpt`` with ``%method ProgExt``).

Covers:
  chemsmart.io.orca.external        .extinp.tmp/.xyz/.engrad files, driver
  ORCAExternalCalculator            launcher, ProgExt path
  ORCAJobSettings.external          ExtOpt routes for opt/TS/IRC/sp
  ORCAInputWriter / ORCAJobRunner   %method block, numerical Hessians
  ``orca ... external`` CLI         reuses opt/ts/... with the backend
  Real ORCA runs                    skipped unless ORCA is installed
"""

import os
import subprocess
import sys
from unittest.mock import MagicMock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.emt import EMT
from ase.units import Bohr, Hartree

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.external import (
    main,
    read_orca_external_input,
    read_xyz_geometry,
    run,
    write_engrad,
)
from chemsmart.jobs.gaussian.external import GaussianExternalCalculator
from chemsmart.jobs.orca.external import ORCAExternalCalculator
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
    ORCANEBJobSettings,
    ORCATSJobSettings,
)
from tests.helpers import (
    make_orca_writer,
    read_engrad,
    run_h2_orca_external,
    run_orca_external_cli,
    write_orca_external_input,
    write_orca_input_text,
)

EMT_CONFIG = {"backend": "ase", "calculator": "ase.calculators.emt.EMT"}

# ── File formats ──────────────────────────────────────────────────────────


class TestORCAExternalFileFormats:
    def test_read_input_resolves_files_against_orca_directory(self, tmp_path):
        path = write_orca_external_input(
            tmp_path,
            ["H", "H"],
            [[0, 0, 0], [0, 0, 0.74]],
            charge=-1,
            multiplicity=2,
            ncores=4,
            gradient=False,
            pointcharges="job_EXT.pc",
        )
        request = read_orca_external_input(path)
        assert request.xyz_file == str(tmp_path / "job_EXT.xyz")
        assert (request.charge, request.multiplicity) == (-1, 2)
        assert request.ncores == 4
        assert request.gradient is False
        assert request.pointcharge_file == str(tmp_path / "job_EXT.pc")
        assert request.engrad_file == str(tmp_path / "job_EXT.engrad")

    def test_point_charge_line_is_optional(self, au2_orca_external_input):
        request = read_orca_external_input(au2_orca_external_input())
        assert request.pointcharge_file is None
        assert request.gradient is True

    def test_engrad_name_without_xyz_extension(self):
        from chemsmart.io.orca.external import ORCAExternalInput

        request = ORCAExternalInput("/run/geom", 0, 1, 1, True)
        assert request.engrad_file == "/run/geom.engrad"

    def test_truncated_input_rejected(self, tmp_path):
        path = tmp_path / "bad.extinp.tmp"
        path.write_text("job_EXT.xyz\n0\n1\n")
        with pytest.raises(ValueError, match="not an ORCA ExtOpt input"):
            read_orca_external_input(str(path))

    def test_xyz_geometry_accepts_labelled_symbols(self, tmp_path):
        xyz = tmp_path / "geom.xyz"
        xyz.write_text("2\ncomment\nC1  0 0 0\ncl 0 0 1.8\n")
        numbers, positions = read_xyz_geometry(str(xyz))
        assert list(numbers) == [6, 17]
        np.testing.assert_allclose(positions[1], [0, 0, 1.8])

    def test_xyz_geometry_with_missing_atoms_rejected(self, tmp_path):
        xyz = tmp_path / "geom.xyz"
        xyz.write_text("3\ncomment\nH 0 0 0\nH 0 0 0.74\n")
        with pytest.raises(ValueError, match="lists 3 atoms"):
            read_xyz_geometry(str(xyz))

    def test_engrad_layout(self, tmp_path):
        path = tmp_path / "job_EXT.engrad"
        gradient = np.arange(6.0).reshape(2, 3) / 100
        write_engrad(path, 2, -1.25, gradient)
        lines = path.read_text().splitlines()
        # natoms, energy and gradient sections, each after 3 comment lines
        assert lines[3].strip() == "2"
        assert float(lines[7]) == pytest.approx(-1.25)
        assert [float(v) for v in lines[11:]] == pytest.approx(
            list(gradient.ravel())
        )
        natoms, energy, parsed = read_engrad(path)
        assert (natoms, energy) == (2, -1.25)
        np.testing.assert_allclose(parsed, gradient)

    def test_energy_only_engrad_keeps_atom_count(self, tmp_path):
        path = tmp_path / "job_EXT.engrad"
        write_engrad(path, 3, -2.0)
        assert read_engrad(path) == (3, -2.0, None)

    def test_engrad_gradient_size_checked(self, tmp_path):
        with pytest.raises(ValueError, match="expected 9"):
            write_engrad(tmp_path / "x.engrad", 3, -1.0, np.zeros((2, 3)))


# ── Driver ────────────────────────────────────────────────────────────────


class TestORCAExternalDriver:
    def test_run_writes_energy_and_gradient_in_atomic_units(
        self, au2_orca_external_input
    ):
        path = au2_orca_external_input()
        message = run(EMT_CONFIG, path)

        reference = Atoms("Au2", positions=[[0, 0, 0], [0, 0, 2.88]])
        reference.calc = EMT()
        natoms, energy, gradient = read_engrad(
            path.replace(".extinp.tmp", ".engrad")
        )
        assert natoms == 2
        # chemsmart uses CODATA 2018 constants, ase.units CODATA 2014
        assert energy == pytest.approx(
            reference.get_potential_energy() / Hartree, rel=1e-7
        )
        np.testing.assert_allclose(
            gradient, -reference.get_forces() * Bohr / Hartree, rtol=1e-7
        )
        assert message.startswith("ase: E = ")

    def test_energy_only_request(self, au2_orca_external_input):
        path = au2_orca_external_input(gradient=False)
        run(EMT_CONFIG, path)
        _, _, gradient = read_engrad(path.replace(".extinp.tmp", ".engrad"))
        assert gradient is None

    def test_point_charges_rejected(self, tmp_path, au2_orca_external_input):
        (tmp_path / "job_EXT.pc").write_text("1\n0.5 0 0 5\n")
        path = au2_orca_external_input(pointcharges="job_EXT.pc")
        with pytest.raises(NotImplementedError, match="Point charges"):
            run(EMT_CONFIG, path)

    def test_empty_point_charge_file_is_ignored(
        self, tmp_path, au2_orca_external_input
    ):
        (tmp_path / "job_EXT.pc").write_text("")
        run(EMT_CONFIG, au2_orca_external_input(pointcharges="job_EXT.pc"))
        assert (tmp_path / "job_EXT.engrad").is_file()

    def test_main_ignores_ext_params(self, au2_orca_external_input, capsys):
        code = main(EMT_CONFIG, [au2_orca_external_input(), "--ignored"])
        assert code == 0
        assert "ase: E = " in capsys.readouterr().out

    def test_main_reports_failure_on_stderr(self, tmp_path, capsys):
        code = main(EMT_CONFIG, [str(tmp_path / "missing.extinp.tmp")])
        assert code == 1
        assert "failed" in capsys.readouterr().err

    def test_main_requires_input_file(self, capsys):
        assert main(EMT_CONFIG, []) == 1
        assert "Usage" in capsys.readouterr().err


# ── Backend configuration ─────────────────────────────────────────────────


class TestORCAExternalCalculator:
    def test_launcher_and_progext(
        self, tmp_path, emt_orca_external_calculator
    ):
        calc = emt_orca_external_calculator
        path = calc.write_script(tmp_path)
        assert os.access(path, os.X_OK)
        assert path == calc.progext(tmp_path)
        assert os.path.isabs(calc.progext("relative/folder"))
        content = open(path).read()
        assert content.startswith(f"#!{sys.executable}\n")
        assert "from chemsmart.io.orca.external import main" in content
        assert "ORCA external launcher" in content

    def test_launcher_name_differs_from_gaussian(self):
        orca_calc = ORCAExternalCalculator(backend="xtb")
        gaussian_calc = GaussianExternalCalculator(backend="xtb")
        assert orca_calc.script_filename.startswith(
            "chemsmart_orca_external_xtb_"
        )
        assert orca_calc.script_filename != gaussian_calc.script_filename
        assert orca_calc != gaussian_calc

    def test_has_no_hessian_options(self, xtb_orca_external_calculator):
        assert xtb_orca_external_calculator.config == {
            "backend": "xtb",
            "calculator": None,
            "options": {"method": "gfn2", "alpb": "water"},
        }
        with pytest.raises(TypeError):
            ORCAExternalCalculator(backend="xtb", hessian="fd")

    def test_launcher_answers_orca_like_orca_does(
        self, tmp_path, emt_orca_external_calculator, au2_orca_external_input
    ):
        """Run the launcher exactly as ORCA would: in ORCA's working
        directory, with the relative .extinp.tmp name and Ext_Params."""
        script = emt_orca_external_calculator.write_script(tmp_path)
        au2_orca_external_input()
        result = subprocess.run(
            [script, "job_EXT.extinp.tmp", "ext-param"],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert result.returncode == 0, result.stderr
        natoms, energy, gradient = read_engrad(tmp_path / "job_EXT.engrad")
        assert natoms == 2
        assert abs(energy) < 10.0  # Hartree, not eV
        assert gradient.shape == (2, 3)


# ── Settings ──────────────────────────────────────────────────────────────


class TestORCAExternalSettings:
    @pytest.mark.parametrize(
        "jobtype, expected",
        [
            ("opt", "! Opt NumFreq ExtOpt"),
            ("ts", "! OptTS NumFreq ExtOpt"),
            ("sp", "! NumFreq ExtOpt"),
        ],
    )
    def test_extopt_replaces_level_of_theory(
        self,
        orca_external_dft_settings,
        xtb_orca_external_calculator,
        jobtype,
        expected,
    ):
        orca_external_dft_settings.jobtype = jobtype
        settings = orca_external_dft_settings.merge(
            {"external": xtb_orca_external_calculator}, keywords=("external",)
        )
        assert settings.route_string == expected

    def test_sp_without_freq(self, xtb_orca_external_calculator):
        settings = ORCAJobSettings(
            jobtype="sp", external=xtb_orca_external_calculator
        )
        assert settings.route_string == "! ExtOpt"

    def test_ts_and_scants_routes(self, xtb_orca_external_calculator):
        settings = ORCATSJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            external=xtb_orca_external_calculator.to_dict(),
        )
        assert settings.external == xtb_orca_external_calculator
        assert settings.route_string == "! OptTS ExtOpt"
        settings.tssearch_type = "scants"
        assert settings.route_string == "! ScanTS ExtOpt"

    def test_irc_route(self, xtb_orca_external_calculator):
        settings = ORCAIRCJobSettings(
            jobtype="irc", external=xtb_orca_external_calculator
        )
        assert settings.route_string.startswith("! IRC")
        assert settings.route_string.endswith("ExtOpt")

    def test_solvation_ignored_and_extra_route_kept(
        self, xtb_orca_external_calculator, caplog
    ):
        settings = ORCAJobSettings(
            jobtype="opt",
            solvent_model="smd",
            solvent_id="water",
            additional_route_parameters="TightOpt ExtOpt",
            external=xtb_orca_external_calculator,
        )
        assert settings.route_string == "! Opt ExtOpt TightOpt"
        assert "solvation is ignored" in caplog.text

    def test_freq_and_numfreq_conflict(self, xtb_orca_external_calculator):
        settings = ORCAJobSettings(
            jobtype="opt",
            freq=True,
            numfreq=True,
            external=xtb_orca_external_calculator,
        )
        with pytest.raises(ValueError, match="both freq and numfreq"):
            settings.route_string

    def test_copy_preserves_external(self, xtb_orca_external_calculator):
        settings = ORCAJobSettings(external=xtb_orca_external_calculator)
        assert settings.copy().external == xtb_orca_external_calculator


# ── Writer and runner ─────────────────────────────────────────────────────


class TestORCAExternalWriter:
    def test_runner_writes_input_and_launcher(
        self,
        tmp_path,
        single_molecule_xyz_file,
        orca_jobrunner_no_scratch,
        emt_orca_external_calculator,
    ):
        settings = ORCAJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            scf_maxiter=500,
            solvent_model="cpcm",
            solvent_id="water",
            external=emt_orca_external_calculator,
        )
        job = ORCAOptJob(
            molecule=Molecule.from_filepath(single_molecule_xyz_file),
            settings=settings,
            label="mol_opt_ase",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)
        orca_jobrunner_no_scratch._prerun(job)
        orca_jobrunner_no_scratch._write_input(job)

        inp = open(orca_jobrunner_no_scratch.job_inputfile).read()
        launcher = emt_orca_external_calculator.progext(tmp_path)
        assert inp.startswith("! Opt ExtOpt\n")
        assert f'%method\n  ProgExt "{launcher}"\nend\n' in inp
        assert "%pal nprocs" in inp
        assert "* xyz 0 1" in inp
        # ORCA's own electronic-structure blocks are not written
        assert "%scf" not in inp
        assert "%cpcm" not in inp
        assert os.access(launcher, os.X_OK)

    def test_ts_hessian_is_numerical(
        self, single_molecule_xyz_file, xtb_orca_external_calculator
    ):
        settings = ORCATSJobSettings(
            charge=0,
            multiplicity=1,
            numhess=False,
            external=xtb_orca_external_calculator,
        )
        writer = make_orca_writer(
            settings,
            molecule=Molecule.from_filepath(single_molecule_xyz_file),
        )
        text = write_orca_input_text(writer)
        assert text.startswith("! OptTS ExtOpt\n")
        assert "Calc_Hess True" in text
        assert "NumHess True" in text
        assert f'ProgExt "{os.path.abspath("/tmp")}/' in text

    @pytest.mark.parametrize(
        "inithess, expected",
        [
            (None, "inithess calc_numfreq"),
            ("calc_anfreq", "inithess calc_numfreq"),
            ("calc_numfreq", "inithess calc_numfreq"),
        ],
    )
    def test_irc_initial_hessian_is_numerical(
        self,
        single_molecule_xyz_file,
        xtb_orca_external_calculator,
        inithess,
        expected,
    ):
        settings = ORCAIRCJobSettings(
            charge=0,
            multiplicity=1,
            inithess=inithess,
            external=xtb_orca_external_calculator,
        )
        writer = make_orca_writer(
            settings,
            molecule=Molecule.from_filepath(single_molecule_xyz_file),
        )
        text = write_orca_input_text(writer)
        assert "%irc\n" in text
        assert expected in text

    def test_neb_rejected(self, xtb_orca_external_calculator):
        settings = ORCANEBJobSettings(
            charge=0, multiplicity=1, external=xtb_orca_external_calculator
        )
        writer = make_orca_writer(settings)
        with pytest.raises(ValueError, match="cannot be combined"):
            write_orca_input_text(writer)


# ── CLI ───────────────────────────────────────────────────────────────────


class TestORCAExternalCLI:
    def test_ts_with_xtb(self, single_molecule_xyz_file):
        job_cls = MagicMock(name="ORCATSJob")
        result, settings = run_orca_external_cli(
            single_molecule_xyz_file,
            "chemsmart.jobs.orca.ts.ORCATSJob",
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
        assert settings.external == ORCAExternalCalculator(
            backend="xtb", options={"method": "gfn2", "alpb": "water"}
        )
        assert settings.route_string.startswith("! OptTS")
        assert settings.route_string.endswith("ExtOpt")
        # ORCA's ts command names OptTS jobs "optts"
        assert job_cls.call_args.kwargs["label"] == "crest_best_optts_xtb"

    def test_opt_with_ase_calculator(self, single_molecule_xyz_file):
        result, settings = run_orca_external_cli(
            single_molecule_xyz_file,
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            "-b",
            "ase",
            "--calculator",
            "ase.calculators.lj.LennardJones",
            "-o",
            "sigma=0.5",
            "--python",
            "/usr/bin/python3",
            "opt",
        )
        assert result.exit_code == 0, result.output
        external = settings.external
        assert external.calculator == "ase.calculators.lj.LennardJones"
        assert external.options == {"sigma": 0.5}
        assert external.python == "/usr/bin/python3"
        assert "ExtOpt" in settings.route_string
        assert "B3LYP" not in settings.route_string.upper()

    def test_user_label_is_kept(self, single_molecule_xyz_file):
        job_cls = MagicMock(name="ORCASinglePointJob")
        result, _ = run_orca_external_cli(
            single_molecule_xyz_file,
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            "-b",
            "xtb",
            "sp",
            orca_args=("-l", "my_label"),
            job_cls=job_cls,
        )
        assert result.exit_code == 0, result.output
        assert job_cls.call_args.kwargs["label"] == "my_label"

    def test_sp_label_has_no_solvent_tag(self, single_molecule_xyz_file):
        """The gas_solv project's sp solvent is ignored by ExtOpt, so it
        must not appear in the label."""
        job_cls = MagicMock(name="ORCASinglePointJob")
        result, settings = run_orca_external_cli(
            single_molecule_xyz_file,
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            "-b",
            "xtb",
            "sp",
            job_cls=job_cls,
        )
        assert result.exit_code == 0, result.output
        assert job_cls.call_args.kwargs["label"] == "crest_best_sp_xtb"
        assert settings.route_string == "! ExtOpt"

    def test_hessian_option_not_offered(self, single_molecule_xyz_file):
        result, _ = run_orca_external_cli(
            single_molecule_xyz_file,
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            "-b",
            "xtb",
            "--hessian",
            "fd",
            "opt",
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "--hessian" in result.output

    def test_ase_backend_requires_calculator(self, single_molecule_xyz_file):
        result, _ = run_orca_external_cli(
            single_molecule_xyz_file,
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            "-b",
            "ase",
            "opt",
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "dotted path" in result.output


# ── Real ORCA runs ────────────────────────────────────────────────────────


class TestORCAExternalRealRun:
    """
    Integration tests that invoke the real ORCA binary.

    Skipped automatically when ORCA is not installed (via the
    ``orca_executable`` fixture). H2 uses ASE's ``LennardJones`` with
    sigma=0.5 Angstrom, as in the Gaussian External real-run tests.
    """

    def test_h2_lj_external_singlepoint(
        self, tmp_path, orca_executable, lj_orca_external_calculator
    ):
        output = run_h2_orca_external(
            orca_executable,
            tmp_path,
            lj_orca_external_calculator,
            "",
            0.74,
            "h2_sp",
        )
        assert "FINAL SINGLE POINT ENERGY" in output

    def test_h2_lj_external_opt(
        self, tmp_path, orca_executable, lj_orca_external_calculator
    ):
        output = run_h2_orca_external(
            orca_executable,
            tmp_path,
            lj_orca_external_calculator,
            "Opt",
            0.90,
            "h2_opt",
        )
        assert "OPTIMIZATION RUN DONE" in output

    def test_h2_lj_external_opt_numfreq(
        self, tmp_path, orca_executable, lj_orca_external_calculator
    ):
        """``NumFreq`` makes ORCA build the Hessian from external
        gradients, as a TS search does for its initial Hessian."""
        output = run_h2_orca_external(
            orca_executable,
            tmp_path,
            lj_orca_external_calculator,
            "Opt NumFreq",
            0.90,
            "h2_numfreq",
        )
        assert "VIBRATIONAL FREQUENCIES" in output
