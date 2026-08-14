from click.testing import CliRunner

from chemsmart.analysis.fukui import discover_fukui_companion_outputs
from chemsmart.cli.fukui import fukui as fukui_analyze
from chemsmart.cli.run import run
from chemsmart.jobs.chain import ChainJob
from chemsmart.jobs.gaussian.fukui import (
    GaussianFukuiJob,
    radical_ion_charge_and_multiplicity,
)
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

        assert isinstance(job, ChainJob)
        assert job.TYPE == "g16fukui"
        assert [phase.name for phase in job.phases] == ["Fukui"]
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

    def test_help_lists_fukui_submit_subcommand(self):
        runner = CliRunner()
        result = runner.invoke(run, ["gaussian", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" in result.output

    def test_help_lists_fukui_analysis_command(self):
        runner = CliRunner()
        result = runner.invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "\n  fukui" in result.output

    def test_submit_help_options(self):
        from chemsmart.cli.gaussian.fukui import fukui as fukui_submit

        runner = CliRunner()
        result = runner.invoke(fukui_submit, ["--help"])
        assert result.exit_code == 0, result.output
        assert "-m, --mode" in result.output
        assert "--radical-cation-charge" in result.output
        assert "--radical-anion-multiplicity" in result.output
        assert "-c, --radical-cation" not in result.output
        assert "-a, --radical-anion" not in result.output

    def test_analyze_help_matches_script_options(self):
        runner = CliRunner()
        result = runner.invoke(fukui_analyze, ["--help"])
        assert result.exit_code == 0, result.output
        assert "-n, --neutral-filename" in result.output
        assert "-c, --radical-cation-filename" in result.output
        assert "-a, --radical-anion-filename" in result.output
        assert "-m, --mode" in result.output
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
