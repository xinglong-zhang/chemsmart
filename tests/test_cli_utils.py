"""Direct unit tests for chemsmart.cli.utils.build_jobs."""

from types import SimpleNamespace
from unittest.mock import MagicMock

from chemsmart.cli.utils import build_jobs


def _make_ctx(molecules, molecule_indices, label, jobrunner=None):
    jobrunner = (
        jobrunner if jobrunner is not None else MagicMock(name="jobrunner")
    )
    return SimpleNamespace(
        obj={
            "jobrunner": jobrunner,
            "molecules": molecules,
            "molecule_indices": molecule_indices,
            "label": label,
        }
    )


class TestBuildJobsSingleMolecule:
    def test_single_molecule_returns_one_job_not_a_list(self):
        molecule = MagicMock(name="molecule")
        ctx = _make_ctx([molecule], None, "mylabel")
        job_cls = MagicMock(name="job_cls")

        result = build_jobs(ctx, job_cls, "settings", True, {})

        job_cls.assert_called_once_with(
            molecule=molecule,
            settings="settings",
            label="mylabel",
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=True,
        )
        assert result is job_cls.return_value

    def test_multiple_molecules_with_no_indices_uses_last_molecule(self):
        mol_a = MagicMock(name="mol_a")
        mol_b = MagicMock(name="mol_b")
        ctx = _make_ctx([mol_a, mol_b], None, "mylabel")
        job_cls = MagicMock(name="job_cls")

        result = build_jobs(ctx, job_cls, "settings", False, {})

        job_cls.assert_called_once_with(
            molecule=mol_b,
            settings="settings",
            label="mylabel",
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=False,
        )
        assert result is job_cls.return_value

    def test_extra_kwargs_are_forwarded_to_job_cls(self):
        molecule = MagicMock(name="molecule")
        ctx = _make_ctx([molecule], None, "mylabel")
        job_cls = MagicMock(name="job_cls")

        build_jobs(ctx, job_cls, "settings", True, {"gfn_version": "gfn2"})

        job_cls.assert_called_once_with(
            molecule=molecule,
            settings="settings",
            label="mylabel",
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=True,
            gfn_version="gfn2",
        )


class TestBuildJobsMultipleMolecules:
    def test_multiple_molecules_with_indices_returns_list_of_jobs(self):
        mol_a = MagicMock(name="mol_a")
        mol_b = MagicMock(name="mol_b")
        ctx = _make_ctx([mol_a, mol_b], [1, 3], "mylabel")
        job_a = MagicMock(name="job_a")
        job_b = MagicMock(name="job_b")
        job_cls = MagicMock(name="job_cls", side_effect=[job_a, job_b])

        result = build_jobs(ctx, job_cls, "settings", True, {})

        assert result == [job_a, job_b]
        assert job_cls.call_count == 2
        job_cls.assert_any_call(
            molecule=mol_a,
            settings="settings",
            label="mylabel_idx1",
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=True,
        )
        job_cls.assert_any_call(
            molecule=mol_b,
            settings="settings",
            label="mylabel_idx3",
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=True,
        )

    def test_multiple_molecules_with_indices_forwards_kwargs_to_each_job(
        self,
    ):
        mol_a = MagicMock(name="mol_a")
        mol_b = MagicMock(name="mol_b")
        ctx = _make_ctx([mol_a, mol_b], [0, 1], "mylabel")
        job_cls = MagicMock(name="job_cls")

        build_jobs(ctx, job_cls, "settings", False, {"nprocs": 4})

        for call in job_cls.call_args_list:
            assert call.kwargs["nprocs"] == 4

    def test_indices_shorter_than_molecules_zips_to_shortest(self):
        # zip() stops at the shorter iterable, so a mismatched
        # molecule_indices list silently caps the number of jobs built
        # rather than raising.
        mol_a = MagicMock(name="mol_a")
        mol_b = MagicMock(name="mol_b")
        ctx = _make_ctx([mol_a, mol_b], [5], "mylabel")
        job_cls = MagicMock(name="job_cls")

        result = build_jobs(ctx, job_cls, "settings", True, {})

        assert len(result) == 1
        job_cls.assert_called_once_with(
            molecule=mol_a,
            settings="settings",
            label="mylabel_idx5",
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=True,
        )
