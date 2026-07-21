"""
Direct unit tests for :class:`GrouperJobRunner`.

Covers directory setup, strategy dispatch (mocking the individual
grouper algorithm classes, which are separately tested via
``test_groupers.py``), success/failure handling in ``_create_process``,
and the XYZ output-writing logic in ``_write_outputs``.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.grouper.job import GrouperJob
from chemsmart.jobs.grouper.runner import GrouperJobRunner


def make_molecules(n=3):
    return [MagicMock(name=f"molecule{i}") for i in range(n)]


@pytest.fixture()
def grouper_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return GrouperJob(
        molecules=make_molecules(3),
        grouping_strategy="rmsd",
        label="mygrp",
        threshold=0.5,
    )


class TestGrouperJobRunnerSetup:
    def test_defaults_scratch_to_false(self, pbs_server):
        runner = GrouperJobRunner(server=pbs_server)
        assert runner.scratch is False

    def test_prerun_creates_output_dir(self, pbs_server, grouper_job):
        runner = GrouperJobRunner(server=pbs_server)
        assert not os.path.exists(grouper_job.output_dir)
        runner._prerun(grouper_job)
        assert os.path.isdir(grouper_job.output_dir)

    def test_prerun_no_op_when_output_dir_exists(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)
        # Should not raise even though the dir already exists.
        runner._prerun(grouper_job)

    def test_get_command_and_executable_are_none(
        self, pbs_server, grouper_job
    ):
        runner = GrouperJobRunner(server=pbs_server)
        assert runner._get_command(grouper_job) is None
        assert runner._get_executable() is None

    def test_run_and_postrun_are_no_ops(self, pbs_server, grouper_job):
        runner = GrouperJobRunner(server=pbs_server)
        runner._run(process=None)
        runner._postrun(grouper_job)


class TestGrouperJobRunnerCreateGrouper:
    @pytest.mark.parametrize(
        "strategy,module,cls_name",
        [
            ("rmsd", "chemsmart.jobs.grouper.rmsd", "BasicRMSDGrouper"),
            ("hrmsd", "chemsmart.jobs.grouper.rmsd", "HungarianRMSDGrouper"),
            ("spyrmsd", "chemsmart.jobs.grouper.rmsd", "SpyRMSDGrouper"),
            ("irmsd", "chemsmart.jobs.grouper.rmsd", "IRMSDGrouper"),
            ("pymolrmsd", "chemsmart.jobs.grouper.rmsd", "PymolRMSDGrouper"),
            (
                "tanimoto",
                "chemsmart.jobs.grouper.tanimoto",
                "TanimotoSimilarityGrouper",
            ),
            (
                "torsion",
                "chemsmart.jobs.grouper.tfd",
                "TorsionFingerprintGrouper",
            ),
            (
                "isomorphism",
                "chemsmart.jobs.grouper.isomorphism",
                "RDKitIsomorphismGrouper",
            ),
            ("formula", "chemsmart.jobs.grouper.formula", "FormulaGrouper"),
            (
                "connectivity",
                "chemsmart.jobs.grouper.connectivity",
                "ConnectivityGrouper",
            ),
            ("energy", "chemsmart.jobs.grouper.energy", "EnergyGrouper"),
        ],
    )
    def test_create_grouper_dispatches_to_correct_class(
        self, pbs_server, grouper_job, strategy, module, cls_name
    ):
        grouper_job.grouping_strategy = strategy
        runner = GrouperJobRunner(server=pbs_server)
        mock_instance = MagicMock()
        with patch(f"{module}.{cls_name}") as mock_cls:
            mock_cls.return_value = mock_instance
            result = runner._create_grouper(grouper_job)

        mock_cls.assert_called_once()
        assert result is mock_instance
        call_kwargs = mock_cls.call_args.kwargs
        assert call_kwargs["molecules"] == grouper_job.molecules
        assert call_kwargs["label"] == grouper_job.label

    def test_create_grouper_unknown_strategy_raises(
        self, pbs_server, grouper_job
    ):
        grouper_job.grouping_strategy = "bogus"
        runner = GrouperJobRunner(server=pbs_server)
        with pytest.raises(ValueError, match="Unknown grouping strategy"):
            runner._create_grouper(grouper_job)

    def test_create_grouper_merges_extra_kwargs(self, pbs_server, grouper_job):
        grouper_job.grouper_kwargs = {"custom_option": "value"}
        runner = GrouperJobRunner(server=pbs_server)
        with patch("chemsmart.jobs.grouper.rmsd.BasicRMSDGrouper") as mock_cls:
            runner._create_grouper(grouper_job)

        assert mock_cls.call_args.kwargs["custom_option"] == "value"


class TestGrouperJobRunnerCreateProcess:
    def test_create_process_success_writes_outputs(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        mock_grouper = MagicMock()
        mock_grouper.group.return_value = ([], [])
        with (
            patch.object(runner, "_create_grouper", return_value=mock_grouper),
            patch.object(runner, "_write_outputs") as mock_write,
        ):
            exit_code = runner._create_process(grouper_job, None, None)

        assert exit_code == 0
        assert grouper_job._grouper is mock_grouper
        mock_write.assert_called_once_with(grouper_job, [], [])

    def test_create_process_failure_writes_errfile(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        with patch.object(
            runner,
            "_create_grouper",
            side_effect=RuntimeError("grouping failed"),
        ):
            exit_code = runner._create_process(grouper_job, None, None)

        assert exit_code == 1
        assert os.path.exists(grouper_job.errfile)
        with open(grouper_job.errfile) as f:
            assert "grouping failed" in f.read()


class SimpleMolecule:
    def __init__(self, energy, symbols, positions):
        self.energy = energy
        self.chemical_symbols = symbols
        self.positions = positions

    @property
    def num_atoms(self):
        return len(self.chemical_symbols)


class TestGrouperJobRunnerWriteOutputs:
    def test_write_outputs_sorts_by_energy_and_writes_xyz(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        mol_high = SimpleMolecule(
            energy=-1.0, symbols=["H"], positions=[(0.0, 0.0, 0.0)]
        )
        mol_low = SimpleMolecule(
            energy=-2.0, symbols=["H"], positions=[(1.0, 1.0, 1.0)]
        )
        groups = [[mol_high, mol_low]]
        group_indices = [[0, 1]]

        runner._write_outputs(grouper_job, groups, group_indices)

        group_file = os.path.join(grouper_job.output_dir, "mygrp_group_1.xyz")
        assert os.path.exists(group_file)
        with open(group_file) as f:
            content = f.read()
        # Lower-energy molecule (mol_low) should be written first.
        assert content.index("Original_Index: 2") < content.index(
            "Original_Index: 1"
        )

    def test_write_outputs_handles_molecules_without_energy(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        runner = GrouperJobRunner(server=pbs_server)

        mol = SimpleMolecule(
            energy=None, symbols=["C"], positions=[(0.0, 0.0, 0.0)]
        )
        groups = [[mol]]
        group_indices = [[0]]

        runner._write_outputs(grouper_job, groups, group_indices)

        group_file = os.path.join(grouper_job.output_dir, "mygrp_group_1.xyz")
        with open(group_file) as f:
            content = f.read()
        assert "Energy: N/A" in content

    def test_write_outputs_uses_conformer_ids_when_provided(
        self, pbs_server, grouper_job
    ):
        os.makedirs(grouper_job.output_dir)
        grouper_job.conformer_ids = ["c1", "c2"]
        runner = GrouperJobRunner(server=pbs_server)

        mol = SimpleMolecule(
            energy=-1.0, symbols=["C"], positions=[(0.0, 0.0, 0.0)]
        )
        groups = [[mol]]
        group_indices = [[1]]

        runner._write_outputs(grouper_job, groups, group_indices)

        group_file = os.path.join(grouper_job.output_dir, "mygrp_group_1.xyz")
        with open(group_file) as f:
            content = f.read()
        assert "Original_Index: c2" in content
