"""
Direct unit tests for the non-multiprocessing logic in
``chemsmart.jobs.iterate.runner``: ``IterateCombination.label``,
``IterateJobRunner._load_molecule``, ``_generate_combinations``,
``_write_outputs``, ``run_single`` (fake mode), and the top-level
``run()`` orchestration.

The multiprocessing-heavy ``run_combinations`` (spawning worker
processes) is intentionally not exercised here — it is integration
behavior better suited to a slower, dedicated test.
"""

import os
from unittest.mock import MagicMock, patch

from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import (
    IterateCombination,
    IterateJobRunner,
    IterateMoleculePool,
)
from chemsmart.jobs.iterate.settings import IterateJobSettings


class TestIterateCombinationLabel:
    def test_label_format(self):
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="benzene",
            skeleton_link_index=5,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="methyl",
            substituent_link_index=1,
        )
        assert combo.label == "benzene_5_methyl_1"


class TestIterateJobRunnerBasics:
    def test_executable_is_none(self):
        runner = IterateJobRunner()
        assert runner.executable is None

    def test_get_command_is_none(self):
        runner = IterateJobRunner()
        assert runner._get_command(job=MagicMock()) is None

    def test_defaults_scratch_to_false(self):
        assert IterateJobRunner.SCRATCH is False
        runner = IterateJobRunner()
        assert runner.scratch is False


class TestRunSingleFakeMode:
    def test_fake_mode_returns_none_without_running(self):
        runner = IterateJobRunner(fake=True)
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        pool = IterateMoleculePool()
        label, result = runner.run_single(combo, pool)
        assert label == combo.label
        assert result is None

    def test_non_fake_mode_delegates_to_task_runner(self):
        runner = IterateJobRunner(fake=False)
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        pool = IterateMoleculePool()
        with patch(
            "chemsmart.jobs.iterate.runner._run_combination_task",
            return_value=("s_1_sub_1", "mol-sentinel"),
        ) as mock_task:
            result = runner.run_single(combo, pool)

        mock_task.assert_called_once_with(combo, pool)
        assert result == ("s_1_sub_1", "mol-sentinel")


class TestLoadMolecule:
    def test_returns_none_without_file_path(self):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {"label": "mymol"}, "skeleton", 0
        )
        assert molecule is None
        assert label == "mymol"

    def test_default_label_when_missing(self):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule({}, "skeleton", 2)
        assert molecule is None
        assert label == "skeleton3"

    def test_loads_molecule_from_file(self, single_molecule_xyz_file):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {"file_path": single_molecule_xyz_file, "label": "myskel"},
            "skeleton",
            0,
        )
        assert molecule is not None
        assert label == "myskel"

    def test_link_index_out_of_bounds_returns_none(
        self, single_molecule_xyz_file
    ):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "link_index": [99999],
            },
            "skeleton",
            0,
        )
        assert molecule is None
        assert label == "myskel"

    def test_skeleton_indices_out_of_bounds_returns_none(
        self, single_molecule_xyz_file
    ):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "skeleton_indices": [99999],
            },
            "skeleton",
            0,
        )
        assert molecule is None
        assert label == "myskel"

    def test_load_failure_returns_none(self):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {"file_path": "/no/such/file.xyz", "label": "bad"},
            "skeleton",
            0,
        )
        assert molecule is None
        assert label == "bad"


class TestGenerateCombinations:
    def test_generates_cross_product_of_valid_configs(
        self, single_molecule_xyz_file
    ):
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "skel1",
                "link_index": [1, 2],
            }
        ]
        settings.substituent_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)

        assert len(pool.skeletons) == 1
        assert len(pool.substituents) == 1
        # 2 skeleton link indices x 1 substituent -> 2 combinations
        assert len(combinations) == 2
        labels = {c.label for c in combinations}
        assert labels == {"skel1_1_sub1_1", "skel1_2_sub1_1"}

    def test_skips_skeleton_without_link_index(self, single_molecule_xyz_file):
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {"file_path": single_molecule_xyz_file, "label": "skel1"}
        ]
        settings.substituent_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert combinations == []

    def test_skips_substituent_without_link_index(
        self, single_molecule_xyz_file
    ):
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "skel1",
                "link_index": [1],
            }
        ]
        settings.substituent_list = [
            {"file_path": single_molecule_xyz_file, "label": "sub1"}
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert combinations == []

    def test_empty_lists_produce_no_combinations(self):
        settings = IterateJobSettings()
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert pool.skeletons == []
        assert pool.substituents == []
        assert combinations == []


class TestWriteOutputs:
    def _make_molecule(self):
        from chemsmart.io.molecules.structure import Molecule

        return Molecule(
            symbols=["H", "H"],
            positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )

    def test_merged_output_writes_all_successful_molecules(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(outputfile="results")
        runner = IterateJobRunner()
        mol = self._make_molecule()
        results = [("combo1", mol), ("combo2", None)]

        runner._write_outputs(results, job)

        assert os.path.exists("results.xyz")
        with open("results.xyz") as f:
            content = f.read()
        assert "combo1" in content
        assert content.count("2\n") == 1  # only one successful molecule

    def test_separate_outputs_writes_one_file_per_result(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(
            separate_outputs=True, output_directory=str(tmp_path / "out")
        )
        runner = IterateJobRunner()
        mol = self._make_molecule()
        results = [("combo1", mol), ("combo2", None)]

        runner._write_outputs(results, job)

        assert os.path.exists(str(tmp_path / "out" / "combo1.xyz"))
        assert not os.path.exists(str(tmp_path / "out" / "combo2.xyz"))

    def test_separate_outputs_falls_back_to_cwd_without_directory(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(separate_outputs=True, output_directory=None)
        runner = IterateJobRunner()
        mol = self._make_molecule()
        results = [("combo1", mol)]

        runner._write_outputs(results, job)

        assert os.path.exists("combo1.xyz")


class TestIterateJobRunnerRun:
    def test_fake_mode_no_op(self):
        runner = IterateJobRunner(fake=True)
        job = MagicMock()
        runner.run(job)
        # No assertions needed beyond "does not raise"; fake mode
        # short-circuits before touching job.settings at all.

    def test_no_valid_combinations_warns_and_returns(self, caplog):
        runner = IterateJobRunner(fake=False)
        job = IterateJob(settings=IterateJobSettings())
        with caplog.at_level("WARNING"):
            runner.run(job)
        assert "No valid combinations" in caplog.text

    def test_full_run_generates_writes_and_uses_run_combinations(
        self, single_molecule_xyz_file, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "skel1",
                "link_index": [1],
            }
        ]
        settings.substituent_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings, outputfile="results")
        runner = IterateJobRunner(fake=False)

        with (
            patch.object(
                runner,
                "run_combinations",
                return_value=[("skel1_1_sub1_1", None)],
            ) as mock_run_combinations,
            patch.object(runner, "_write_outputs") as mock_write,
        ):
            runner.run(job)

        mock_run_combinations.assert_called_once()
        mock_write.assert_called_once()
