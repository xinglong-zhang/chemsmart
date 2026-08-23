"""
Direct unit tests for :class:`GrouperJob`.

The grouping *algorithms* themselves (RMSD, Tanimoto, TFD, etc.) are
already thoroughly tested in ``test_groupers.py``. These tests cover
the ``GrouperJob`` wrapper class: constructor validation, path
properties, completion checks, and ``from_filename``.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.grouper.job import GrouperJob


def make_molecules(n=3):
    return [MagicMock(name=f"molecule{i}") for i in range(n)]


class TestGrouperJobConstruction:
    def test_rejects_non_list_molecules(self):
        with pytest.raises(ValueError, match="must be a list"):
            GrouperJob(molecules="not-a-list", grouping_strategy="rmsd")

    def test_rejects_fewer_than_two_molecules(self):
        with pytest.raises(ValueError, match="at least 2"):
            GrouperJob(molecules=make_molecules(1), grouping_strategy="rmsd")

    def test_label_defaults_to_grouper(self):
        job = GrouperJob(molecules=make_molecules(2), grouping_strategy="rmsd")
        assert job.label == "grouper"

    def test_custom_label_used(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="my_group",
        )
        assert job.label == "my_group"

    def test_stores_configuration(self):
        job = GrouperJob(
            molecules=make_molecules(3),
            grouping_strategy="tanimoto",
            threshold=0.9,
            num_groups=2,
            ignore_hydrogens=True,
            num_procs=4,
            matrix_format="csv",
            energy_type="G",
        )
        assert job.grouping_strategy == "tanimoto"
        assert job.threshold == 0.9
        assert job.num_groups == 2
        assert job.ignore_hydrogens is True
        assert job.num_procs == 4
        assert job.matrix_format == "csv"
        assert job.energy_type == "G"
        assert job.num_molecules == 3


class TestGrouperJobPaths:
    def test_job_basename_returns_label(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.job_basename == "mygrp"

    def test_output_dir_uses_label(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.output_dir == os.path.join(job.folder, "mygrp_group_result")

    def test_output_dir_without_label_falls_back(self):
        job = GrouperJob(molecules=make_molecules(2), grouping_strategy="rmsd")
        job.label = None
        assert job.output_dir == os.path.join(job.folder, "group_result")

    def test_outputfile_path(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.outputfile == os.path.join(
            job.output_dir, "mygrp_group_1.xyz"
        )

    def test_logfile_and_errfile_paths(self):
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        assert job.logfile == os.path.join(job.folder, "log.mygrp")
        assert job.errfile == os.path.join(job.folder, "mygrp.err")


class TestGrouperJobCompletion:
    def test_not_complete_when_output_dir_missing(self):
        job = GrouperJob(molecules=make_molecules(2), grouping_strategy="rmsd")
        assert job.is_complete() is False

    def test_complete_when_output_file_exists(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            label="mygrp",
        )
        os.makedirs(job.output_dir)
        with open(job.outputfile, "w") as f:
            f.write("1\ncomment\nC 0 0 0\n")
        assert job.is_complete() is True


class TestGrouperJobRun:
    def test_run_delegates_to_jobrunner(self):
        mock_runner = MagicMock()
        job = GrouperJob(
            molecules=make_molecules(2),
            grouping_strategy="rmsd",
            jobrunner=mock_runner,
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestGrouperJobFromFilename:
    def test_from_filename_rejects_too_few_molecules(
        self, single_molecule_xyz_file
    ):
        with patch(
            "chemsmart.jobs.grouper.job.Molecule.from_filepath",
            return_value=[MagicMock()],
        ):
            with pytest.raises(ValueError, match="at least 2"):
                GrouperJob.from_filename(
                    filename=single_molecule_xyz_file,
                    grouping_strategy="rmsd",
                )

    def test_from_filename_builds_job_with_jobrunner(
        self, multiple_molecules_xyz_file
    ):
        mock_jobrunner = MagicMock()
        job = GrouperJob.from_filename(
            filename=multiple_molecules_xyz_file,
            grouping_strategy="rmsd",
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, GrouperJob)
        assert job.jobrunner is mock_jobrunner
        assert job.num_molecules >= 2

    def test_from_filename_derives_label_from_filename(
        self, multiple_molecules_xyz_file
    ):
        mock_jobrunner = MagicMock()
        job = GrouperJob.from_filename(
            filename=multiple_molecules_xyz_file,
            grouping_strategy="rmsd",
            jobrunner=mock_jobrunner,
        )
        expected_label = os.path.basename(multiple_molecules_xyz_file).split(
            "."
        )[0]
        assert job.label == expected_label

    def test_from_filename_creates_jobrunner_when_not_given(
        self, multiple_molecules_xyz_file
    ):
        mock_runner_instance = MagicMock()
        with patch(
            "chemsmart.jobs.runner.JobRunner.from_job",
            return_value=mock_runner_instance,
        ) as mock_from_job:
            job = GrouperJob.from_filename(
                filename=multiple_molecules_xyz_file,
                grouping_strategy="rmsd",
            )

        mock_from_job.assert_called_once()
        assert job.jobrunner is mock_runner_instance
