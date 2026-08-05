"""
Direct unit tests for :class:`PyMOLIRCMovieJob` (``chemsmart.jobs.mol.irc``).

CLI-level invocation (through ``chemsmart run mol irc``) is covered by
``test_mol_cli.py``, but only patches the class entirely, so
``from_files``'s branching logic (reactant+product vs. full-trajectory
vs. invalid combinations) was never actually exercised.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.mol.irc import PyMOLIRCMovieJob


class TestPyMOLIRCMovieJobFromFiles:
    def test_no_files_provided_raises(self):
        with pytest.raises(ValueError, match="No reactant or product"):
            PyMOLIRCMovieJob.from_files(
                reactant_file=None,
                product_file=None,
                all_file=None,
                label="x",
            )

    def test_all_three_files_provided_raises(self):
        with pytest.raises(ValueError, match="only provide reactant"):
            PyMOLIRCMovieJob.from_files(
                reactant_file="r.log",
                product_file="p.log",
                all_file="full.log",
                label="x",
            )

    def test_only_reactant_file_raises(self):
        with pytest.raises(ValueError, match="only provide reactant"):
            PyMOLIRCMovieJob.from_files(
                reactant_file="r.log",
                product_file=None,
                all_file=None,
                label="x",
            )

    def test_reactant_and_product_files_build_job(self):
        reactant_mols = [MagicMock(name="r1"), MagicMock(name="r2")]
        product_mols = [MagicMock(name="p1"), MagicMock(name="p2")]
        mock_jobrunner = MagicMock()

        with patch(
            "chemsmart.jobs.mol.irc.Molecule.from_filepath"
        ) as mock_from_filepath:
            mock_from_filepath.side_effect = [reactant_mols, product_mols]
            job = PyMOLIRCMovieJob.from_files(
                reactant_file="path/to/reactant_ircr.log",
                product_file="path/to/reactant_ircf.log",
                all_file=None,
                label="ignored_label",
                jobrunner=mock_jobrunner,
            )

        assert isinstance(job, PyMOLIRCMovieJob)
        assert job.jobrunner is mock_jobrunner
        assert job.molecule == reactant_mols + product_mols
        # label derived from common prefix of the two basenames, not the
        # explicitly passed "ignored_label"
        assert job.label.startswith("reactant_irc")

    def test_reactant_and_product_files_read_in_correct_order_and_index(
        self,
    ):
        with patch(
            "chemsmart.jobs.mol.irc.Molecule.from_filepath"
        ) as mock_from_filepath:
            mock_from_filepath.side_effect = [[], []]
            PyMOLIRCMovieJob.from_files(
                reactant_file="reactant.log",
                product_file="product.log",
                all_file=None,
                label="x",
                jobrunner=MagicMock(),
            )

        assert mock_from_filepath.call_args_list[0].args[0] == "reactant.log"
        assert mock_from_filepath.call_args_list[0].kwargs["index"] == "::-1"
        assert mock_from_filepath.call_args_list[1].args[0] == "product.log"
        assert mock_from_filepath.call_args_list[1].kwargs["index"] == ":"

    def test_all_file_builds_job(self):
        all_mols = [MagicMock(name="m1"), MagicMock(name="m2")]
        mock_jobrunner = MagicMock()

        with patch(
            "chemsmart.jobs.mol.irc.Molecule.from_filepath",
            return_value=all_mols,
        ) as mock_from_filepath:
            job = PyMOLIRCMovieJob.from_files(
                reactant_file=None,
                product_file=None,
                all_file="path/to/full_irc.log",
                label="ignored_label",
                jobrunner=mock_jobrunner,
            )

        mock_from_filepath.assert_called_once_with(
            "path/to/full_irc.log", index=":"
        )
        assert job.molecule == all_mols
        assert job.label == "full_irc_movie"

    def test_from_files_creates_jobrunner_when_not_given(self):
        mock_runner_instance = MagicMock()
        with (
            patch(
                "chemsmart.jobs.mol.irc.Molecule.from_filepath",
                return_value=[],
            ),
            patch(
                "chemsmart.jobs.runner.JobRunner.from_job",
                return_value=mock_runner_instance,
            ) as mock_from_job,
        ):
            job = PyMOLIRCMovieJob.from_files(
                reactant_file=None,
                product_file=None,
                all_file="full.log",
                label="x",
            )

        mock_from_job.assert_called_once()
        assert job.jobrunner is mock_runner_instance


class TestPyMOLIRCMovieJobConstruction:
    def test_label_suffixed_with_movie(self):
        job = PyMOLIRCMovieJob(molecules=[MagicMock()], label="irc_run")
        assert job.label == "irc_run_movie"

    def test_type_identifier(self):
        job = PyMOLIRCMovieJob(molecules=[MagicMock()], label="irc_run")
        assert job.TYPE == "pymol_ircmovie"
