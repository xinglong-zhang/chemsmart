from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    build_molecule,
    build_orca_settings,
    dry_run_input,
)
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
)


class TestBuildSettingsAndJob:
    def test_build_gaussian_settings_returns_settings(self):
        settings = build_gaussian_settings("B3LYP", "6-31G*")

        assert isinstance(settings, GaussianJobSettings)
        assert settings.functional == "B3LYP"
        assert settings.basis == "6-31G*"

    def test_build_gaussian_settings_exposes_frequency_flags(self):
        settings = build_gaussian_settings(
            "B3LYP",
            "6-31G*",
            freq=True,
            numfreq=False,
            additional_opt_options_in_route="maxstep=8",
            additional_route_parameters="scf=tight",
        )

        assert settings.freq is True
        assert settings.numfreq is False
        assert settings.additional_opt_options_in_route == "maxstep=8"
        assert settings.additional_route_parameters == "scf=tight"

    def test_build_gaussian_opt_job(
        self,
        single_molecule_xyz_file,
    ):
        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings("B3LYP", "6-31G*")

        job = build_job("gaussian.opt", molecule=molecule, settings=settings)

        assert isinstance(job, GaussianOptJob)
        assert job.molecule.empirical_formula == molecule.empirical_formula
        assert job.settings.functional == "B3LYP"
        assert job.settings.basis == "6-31G*"
        assert job.settings.jobtype == "opt"

    def test_unknown_job_kind_lists_supported_kinds(
        self,
        single_molecule_xyz_file,
    ):
        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings("B3LYP", "6-31G*")

        with pytest.raises(ValueError) as excinfo:
            build_job("gaussian.unknown", molecule=molecule, settings=settings)

        message = str(excinfo.value)
        assert "gaussian.opt" in message
        assert "gaussian.ts" in message
        assert "gaussian.freq" in message
        assert "orca.opt" in message

    def test_build_orca_opt_job(
        self,
        single_molecule_xyz_file,
    ):
        try:
            from chemsmart.jobs.orca.opt import ORCAOptJob
        except Exception as exc:
            pytest.skip(f"ORCA jobs are not importable: {exc}")

        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_orca_settings("B3LYP", "def2-SVP")

        job = build_job("orca.opt", molecule=molecule, settings=settings)

        assert isinstance(job, ORCAOptJob)
        assert job.settings.functional == "B3LYP"
        assert job.settings.basis == "def2-SVP"
        assert job.settings.jobtype == "opt"

    def test_build_gaussian_irc_job_promotes_irc_settings(
        self,
        single_molecule_xyz_file,
    ):
        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings("B3LYP", "6-31G*")

        job = build_job("gaussian.irc", molecule=molecule, settings=settings)

        assert isinstance(job, GaussianIRCJob)
        assert isinstance(job.settings, GaussianIRCJobSettings)
        assert job.settings.jobtype == "irc"
        assert "irc=(" in job.settings.route_string
        assert "stepsize=10" in job.settings.route_string
        assert "maxpoints=20" in job.settings.route_string

    def test_build_gaussian_opt_job_with_freq_writes_single_route(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
    ):
        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings(
            "B3LYP",
            "6-31G*",
            freq=True,
        )

        job = build_job(
            "gaussian.opt",
            molecule=molecule,
            settings=settings,
            label="gaussian_opt_freq",
        )
        job.set_folder(str(tmp_path))

        result = dry_run_input(job)

        route_line = next(
            line.strip()
            for line in result["content"].splitlines()
            if line.lstrip().startswith("#")
        ).lower()
        assert " opt" in route_line
        assert " freq" in route_line
        assert Path(result["inputfile"]).name == "gaussian_opt_freq.com"
        assert sorted(path.name for path in tmp_path.glob("*.com")) == [
            "gaussian_opt_freq.com"
        ]

    def test_scan_with_definition_produces_modredundant_section(
        self,
        tmp_path: Path,
        single_molecule_xyz_file,
    ):
        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings(
            "B3LYP",
            "6-31G*",
            scan_definition="B 1 2 S 10 0.05",
        )

        job = build_job(
            "gaussian.scan",
            molecule=molecule,
            settings=settings,
            label="gaussian_scan",
        )
        job.set_folder(str(tmp_path))

        result = dry_run_input(job)

        assert "# opt=modredundant" in result["content"].lower()
        assert "B 1 2 S 10 0.05" in result["content"]

    def test_scan_without_definition_raises_value_error(
        self,
        single_molecule_xyz_file,
    ):
        molecule = build_molecule(single_molecule_xyz_file)
        settings = build_gaussian_settings("B3LYP", "6-31G*")

        with pytest.raises(ValueError) as excinfo:
            build_job("gaussian.scan", molecule=molecule, settings=settings)

        assert "scan_definition" in str(excinfo.value)
