"""Job factory methods that select their own job runner.

``server``, ``scratch`` and ``fake`` once reached ``JobRunner.from_job`` twice,
so a factory raised ``TypeError`` whenever one of them was supplied.
"""

from unittest.mock import patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import (
    GaussianComJob,
    GaussianGeneralJob,
    GaussianJob,
)
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.grouper.job import GrouperJob
from chemsmart.jobs.mol.irc import PyMOLIRCMovieJob
from chemsmart.jobs.mol.visualize import PyMOLVisualizationJob
from chemsmart.jobs.orca.job import ORCAGeneralJob, ORCAInpJob, ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.thermochemistry.boltzmann import (
    BoltzmannAverageThermochemistryJob,
)
from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
from tests.helpers import assert_runner_was_configured


class TestGaussianFactories:
    def test_from_filename(
        self, single_molecule_xyz_file, runner_options, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianOptJob.from_filename(
            single_molecule_xyz_file,
            settings=GaussianJobSettings(),
            label="from_file",
            **runner_options,
        )
        assert_runner_was_configured(job, runner_options)

    def test_from_pubchem(self, single_molecule_xyz_file, runner_options):
        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        with patch.object(Molecule, "from_pubchem", return_value=molecule):
            job = GaussianOptJob.from_pubchem(
                "962",
                settings=GaussianJobSettings(),
                label="fetched",
                **runner_options,
            )
        assert_runner_was_configured(job, runner_options)

    @pytest.mark.xfail(
        strict=True,
        raises=ValueError,
        reason="from_filename passes molecule=None, which the job rejects",
    )
    def test_com_job_from_filename(
        self, gaussian_frozen_opt_inputfile, runner_options
    ):
        job = GaussianComJob.from_filename(
            gaussian_frozen_opt_inputfile, label="com", **runner_options
        )
        assert_runner_was_configured(job, runner_options)

    @pytest.mark.parametrize(
        ("jobtype", "job_class"),
        [
            ("opt", GaussianOptJob),
            ("com", GaussianComJob),
            ("g16", GaussianGeneralJob),
        ],
    )
    def test_from_jobtype(
        self, jobtype, job_class, single_molecule_xyz_file, runner_options
    ):
        job = GaussianJob.from_jobtype(
            jobtype,
            Molecule.from_filepath(single_molecule_xyz_file),
            settings=GaussianJobSettings(),
            label="typed",
            **runner_options,
        )
        assert isinstance(job, job_class)
        assert_runner_was_configured(job, runner_options)


class TestORCAFactories:
    def test_from_filename(self, single_molecule_xyz_file, runner_options):
        job = ORCAOptJob.from_filename(
            single_molecule_xyz_file,
            settings=ORCAJobSettings(),
            label="from_file",
            **runner_options,
        )
        assert_runner_was_configured(job, runner_options)

    def test_from_pubchem(self, single_molecule_xyz_file, runner_options):
        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        with patch.object(Molecule, "from_pubchem", return_value=molecule):
            job = ORCAOptJob.from_pubchem(
                "962",
                settings=ORCAJobSettings(),
                label="fetched",
                **runner_options,
            )
        assert_runner_was_configured(job, runner_options)

    def test_inp_job_from_filename(self, water_sp_input_path, runner_options):
        job = ORCAInpJob.from_filename(
            water_sp_input_path, label="inp", **runner_options
        )
        assert_runner_was_configured(job, runner_options)

    @pytest.mark.parametrize(
        ("jobtype", "job_class"),
        [
            ("opt", ORCAOptJob),
            ("inp", ORCAInpJob),
            ("orca", ORCAGeneralJob),
        ],
    )
    def test_from_jobtype(
        self, jobtype, job_class, single_molecule_xyz_file, runner_options
    ):
        job = ORCAJob.from_jobtype(
            jobtype,
            Molecule.from_filepath(single_molecule_xyz_file),
            settings=ORCAJobSettings(),
            label="typed",
            **runner_options,
        )
        assert isinstance(job, job_class)
        assert_runner_was_configured(job, runner_options)


class TestPyMOLFactories:
    def test_from_filename(self, single_molecule_xyz_file, runner_options):
        job = PyMOLVisualizationJob.from_filename(
            single_molecule_xyz_file, label="view", **runner_options
        )
        assert_runner_was_configured(job, runner_options)

    def test_from_pubchem(self, single_molecule_xyz_file, runner_options):
        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        with patch.object(Molecule, "from_pubchem", return_value=molecule):
            job = PyMOLVisualizationJob.from_pubchem(
                "962", label="fetched", **runner_options
            )
        assert_runner_was_configured(job, runner_options)

    def test_irc_movie_from_files(
        self, multiple_molecules_xyz_file, runner_options
    ):
        job = PyMOLIRCMovieJob.from_files(
            reactant_file=None,
            product_file=None,
            all_file=multiple_molecules_xyz_file,
            label="irc",
            **runner_options,
        )
        assert_runner_was_configured(job, runner_options)


class TestOtherFactories:
    def test_grouper_from_filename(
        self, multiple_molecules_xyz_file, runner_options
    ):
        job = GrouperJob.from_filename(
            multiple_molecules_xyz_file,
            grouping_strategy="formula",
            index=":",
            **runner_options,
        )
        assert_runner_was_configured(job, runner_options)

    def test_thermochemistry_from_filename(
        self, gaussian_co2_opt_outfile, runner_options
    ):
        job = ThermochemistryJob.from_filename(
            gaussian_co2_opt_outfile, **runner_options
        )
        assert_runner_was_configured(job, runner_options)

    def test_boltzmann_from_files(
        self, gaussian_co2_opt_outfile, runner_options
    ):
        job = BoltzmannAverageThermochemistryJob.from_files(
            files=[gaussian_co2_opt_outfile],
            filename=gaussian_co2_opt_outfile,
            **runner_options,
        )
        assert_runner_was_configured(job, runner_options)


def test_extra_job_arguments_reach_the_job(
    single_molecule_xyz_file, runner_options
):
    job = GaussianOptJob.from_filename(
        single_molecule_xyz_file,
        settings=GaussianJobSettings(),
        label="extra",
        skip_completed=False,
        **runner_options,
    )
    assert job.skip_completed is False


def test_a_supplied_runner_is_used_as_is(
    single_molecule_xyz_file, gaussian_jobrunner_no_scratch, pbs_server
):
    job = GaussianOptJob.from_filename(
        single_molecule_xyz_file,
        settings=GaussianJobSettings(),
        label="given",
        jobrunner=gaussian_jobrunner_no_scratch,
        server=pbs_server,
    )
    assert job.jobrunner is gaussian_jobrunner_no_scratch
