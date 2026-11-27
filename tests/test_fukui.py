import pytest
from click.testing import CliRunner

from chemsmart.analysis.fukui import (
    analyze_fukui,
    discover_fukui_companion_outputs,
    radical_ion_charge_and_multiplicity,
)
from chemsmart.cli.fukui import fukui as fukui_analyze
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.gaussian.fukui import GaussianFukuiJob
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.runner import FakeGaussianJobRunner
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.wbi import GaussianWBIJob
from chemsmart.jobs.runner import JobRunner


class TestFukuiHelpers:
    def test_closed_shell_ions_are_doublets(self):
        assert radical_ion_charge_and_multiplicity(0, 1, +1) == (1, 2)
        assert radical_ion_charge_and_multiplicity(0, 1, -1) == (-1, 2)

    def test_open_shell_ions_decrease_multiplicity(self):
        assert radical_ion_charge_and_multiplicity(0, 2, +1) == (1, 1)

    def test_discover_companions_from_neutral_label(self, tmp_path):
        neutral = tmp_path / "mol_n.log"
        cation = tmp_path / "mol_rc.log"
        anion = tmp_path / "mol_ra.out"
        neutral.write_text("n")
        cation.write_text("c")
        anion.write_text("a")

        found = discover_fukui_companion_outputs(str(neutral))
        assert found["radical_cation"] == str(cation)
        assert found["radical_anion"] == str(anion)


class DummyFukuiOutput:
    def __init__(self, energy, charges):
        self.energies = [energy]
        self.mulliken_atomic_charges = charges


class TestAnalyzeFukuiOutput:
    def test_writes_results_table(self, tmp_path, monkeypatch, caplog):
        import logging

        caplog.set_level(logging.INFO, logger="chemsmart.analysis.fukui")
        outputs = {
            "n.log": DummyFukuiOutput(-100.0, {"C1": 0.10, "H2": 0.00}),
            "c.log": DummyFukuiOutput(-99.5, {"C1": 0.40, "H2": 0.10}),
            "a.log": DummyFukuiOutput(-100.4, {"C1": -0.20, "H2": -0.10}),
        }

        def fake_load(filename):
            return outputs[filename], "gaussian"

        monkeypatch.setattr("chemsmart.analysis.fukui._load_output", fake_load)
        output_path = tmp_path / "subdir" / "fukui.dat"
        analyze_fukui(
            neutral_filename="n.log",
            radical_cation_filename="c.log",
            radical_anion_filename="a.log",
            mode="mulliken",
            output=str(output_path),
        )
        text = output_path.read_text()
        assert "Ionization energy = 0.5" in text
        assert "Fukui Minus (f-)" in text
        assert "C1" in text
        assert "Ionization energy" not in caplog.text

    def test_nbo_rejects_orca_outputs(self, monkeypatch):
        def fake_load(filename):
            return DummyFukuiOutput(-100.0, {"C1": 0.10}), "orca"

        monkeypatch.setattr("chemsmart.analysis.fukui._load_output", fake_load)
        with pytest.raises(
            ValueError, match="NBO charges are only available for Gaussian"
        ):
            analyze_fukui(
                neutral_filename="n.out",
                radical_cation_filename="c.out",
                radical_anion_filename="a.out",
                mode="nbo",
            )

    def test_cli_passes_output(self, tmp_path, mocker, caplog):
        import logging

        caplog.set_level(logging.INFO, logger="chemsmart.cli.fukui")
        mock = mocker.patch("chemsmart.cli.fukui.analyze_fukui")
        neutral = tmp_path / "mol_n.log"
        cation = tmp_path / "mol_rc.log"
        anion = tmp_path / "mol_ra.log"
        neutral.write_text("n")
        cation.write_text("c")
        anion.write_text("a")
        output_path = tmp_path / "fukui.dat"
        runner = CliRunner()
        result = runner.invoke(
            fukui_analyze,
            ["-n", str(neutral), "-o", str(output_path)],
        )
        assert result.exit_code == 0, result.output
        mock.assert_called_once()
        assert mock.call_args.kwargs["output"] == str(output_path)
        assert mock.call_args.kwargs["radical_cation_filename"] == str(cation)
        assert mock.call_args.kwargs["radical_anion_filename"] == str(anion)
        assert "Auto-discovered" not in caplog.text


class TestGaussianFukuiJob:
    def test_creates_neutral_cation_and_anion_jobs(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        settings = GaussianJobSettings(
            functional="B3LYP",
            basis="6-31G*",
            charge=0,
            multiplicity=1,
        )
        job = GaussianFukuiJob(
            molecule=mol,
            settings=settings,
            label="phenol_fukui",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianJob)
        assert job.TYPE == "g16fukui"
        assert [phase.name for phase in job.phases] == ["Fukui"]
        fukui_phase = job.phase_by_name("Fukui")
        assert fukui_phase.stop_on_incomplete is True
        assert fukui_phase.require_complete is True
        assert isinstance(job.neutral_job, GaussianSinglePointJob)
        assert job.neutral_job.label == "phenol_fukui_n"
        assert job.cation_job.label == "phenol_fukui_rc"
        assert job.anion_job.label == "phenol_fukui_ra"
        assert job.cation_job.settings.charge == 1
        assert job.cation_job.settings.multiplicity == 2
        assert job.anion_job.settings.charge == -1
        assert job.anion_job.settings.multiplicity == 2

    def test_nbo_mode_uses_wbi_jobs(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        settings = GaussianJobSettings(
            functional="B3LYP",
            basis="6-31G*",
            charge=0,
            multiplicity=1,
        )
        job = GaussianFukuiJob(
            molecule=mol,
            settings=settings,
            label="phenol_fukui",
            jobrunner=gaussian_jobrunner_no_scratch,
            mode="nbo",
        )

        assert isinstance(job.neutral_job, GaussianWBIJob)
        assert job.neutral_job.settings.jobtype == "wbi"

    def test_user_overrides_ion_charge_and_multiplicity(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = GaussianJobSettings(
            functional="B3LYP",
            basis="6-31G*",
            charge=0,
            multiplicity=1,
        )
        job = GaussianFukuiJob(
            molecule=mol,
            settings=settings,
            label="phenol_fukui",
            jobrunner=gaussian_jobrunner_no_scratch,
            radical_cation_charge=2,
            radical_cation_multiplicity=1,
            radical_anion_charge=-2,
            radical_anion_multiplicity=1,
        )

        assert job.cation_job.settings.charge == 2
        assert job.cation_job.settings.multiplicity == 1
        assert job.anion_job.settings.charge == -2
        assert job.anion_job.settings.multiplicity == 1


class TestORCAFukuiJob:
    def test_creates_neutral_cation_and_anion_jobs(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.fukui import ORCAFukuiJob
        from chemsmart.jobs.orca.job import ORCAJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings
        from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        job = ORCAFukuiJob(
            molecule=mol,
            settings=settings,
            label="phenol_fukui",
            jobrunner=orca_jobrunner_no_scratch,
        )

        assert isinstance(job, ORCAJob)
        assert job.TYPE == "orcafukui"
        assert [phase.name for phase in job.phases] == ["Fukui"]
        fukui_phase = job.phase_by_name("Fukui")
        assert fukui_phase.stop_on_incomplete is True
        assert fukui_phase.require_complete is True
        assert isinstance(job.neutral_job, ORCASinglePointJob)
        assert job.neutral_job.label == "phenol_fukui_n"
        assert job.cation_job.label == "phenol_fukui_rc"
        assert job.anion_job.label == "phenol_fukui_ra"
        assert job.cation_job.settings.charge == 1
        assert job.cation_job.settings.multiplicity == 2
        assert job.anion_job.settings.charge == -1
        assert job.anion_job.settings.multiplicity == 2

    def test_hirshfeld_mode_adds_route_keyword(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.fukui import ORCAFukuiJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        job = ORCAFukuiJob(
            molecule=mol,
            settings=settings,
            label="phenol_fukui",
            jobrunner=orca_jobrunner_no_scratch,
            mode="hirshfeld",
        )

        assert (
            job.neutral_job.settings.additional_route_parameters == "Hirshfeld"
        )
        assert "Hirshfeld" in job.neutral_job.settings.route_string

    def test_nbo_and_cm5_modes_are_rejected(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.fukui import ORCAFukuiJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        for mode in ("nbo", "cm5"):
            with pytest.raises(ValueError, match="Supported ORCA modes"):
                ORCAFukuiJob(
                    molecule=mol,
                    settings=settings,
                    label="phenol_fukui",
                    jobrunner=orca_jobrunner_no_scratch,
                    mode=mode,
                )

    def test_user_overrides_ion_charge_and_multiplicity(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.fukui import ORCAFukuiJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        settings = ORCAJobSettings(
            functional="B3LYP",
            basis="def2-SVP",
            charge=0,
            multiplicity=1,
        )
        job = ORCAFukuiJob(
            molecule=mol,
            settings=settings,
            label="phenol_fukui",
            jobrunner=orca_jobrunner_no_scratch,
            radical_cation_charge=2,
            radical_cation_multiplicity=1,
            radical_anion_charge=-2,
            radical_anion_multiplicity=1,
        )

        assert job.cation_job.settings.charge == 2
        assert job.cation_job.settings.multiplicity == 1
        assert job.anion_job.settings.charge == -2
        assert job.anion_job.settings.multiplicity == 1


class TestFukuiCLI:
    def test_gaussian_runner_accepts_fukui_type(self, pbs_server):
        from types import SimpleNamespace

        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="g16fukui"),
            server=pbs_server,
            scratch=False,
            fake=True,
        )
        assert isinstance(runner, FakeGaussianJobRunner)

    def test_orca_runner_accepts_fukui_type(self, pbs_server):
        from types import SimpleNamespace

        from chemsmart.jobs.orca.runner import FakeORCAJobRunner

        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcafukui"),
            server=pbs_server,
            scratch=False,
            fake=True,
        )
        assert isinstance(runner, FakeORCAJobRunner)

    def test_help_lists_fukui_submit_subcommand(self):
        runner = CliRunner()
        result = runner.invoke(run, ["gaussian", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" in result.output
        result = runner.invoke(run, ["orca", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" in result.output

    def test_help_lists_fukui_analysis_command(self):
        runner = CliRunner()
        result = runner.invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" in result.output
        result = runner.invoke(run, ["chain", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" in result.output

    def test_sub_help_does_not_list_analysis_fukui(self):
        runner = CliRunner()
        result = runner.invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" not in result.output

    def test_submit_help_via_run(self, single_molecule_xyz_file):
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--fake",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "fukui",
                "--help",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "--mode" in result.output
        assert "-m, --mode" not in result.output

    def test_submit_help_options(self):
        from chemsmart.cli.gaussian.fukui import fukui as gaussian_fukui_submit
        from chemsmart.cli.orca.fukui import fukui as orca_fukui_submit

        runner = CliRunner()
        result = runner.invoke(gaussian_fukui_submit, ["--help"])
        assert result.exit_code == 0, result.output
        assert "--mode" in result.output
        assert "-m, --mode" not in result.output
        assert "--radical-cation-charge" in result.output
        assert "--radical-anion-multiplicity" in result.output
        assert "-c, --radical-cation" not in result.output
        assert "-a, --radical-anion" not in result.output

        result = runner.invoke(orca_fukui_submit, ["--help"])
        assert result.exit_code == 0, result.output
        assert "--mode" in result.output
        assert "-m, --mode" not in result.output
        assert "--radical-cation-charge" in result.output
        assert "--radical-anion-multiplicity" in result.output

    def test_analyze_help_matches_script_options(self):
        runner = CliRunner()
        result = runner.invoke(fukui_analyze, ["--help"])
        assert result.exit_code == 0, result.output
        assert "-n, --neutral-filename" in result.output
        assert "-c, --radical-cation-filename" in result.output
        assert "-a, --radical-anion-filename" in result.output
        assert "-m, --mode" in result.output
        assert "-o, --output" in result.output
        assert "Charges to be used for Fukui Indices" in result.output

    def test_analyze_requires_ion_file(self, tmp_path):
        neutral = tmp_path / "only_n.log"
        neutral.write_text("n")
        runner = CliRunner()
        result = runner.invoke(
            fukui_analyze, ["-n", str(neutral), "-m", "mulliken"]
        )
        assert result.exit_code != 0
        assert (
            "radical-cation" in result.output
            or "radical-anion" in result.output
        )
