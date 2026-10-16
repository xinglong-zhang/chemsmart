"""
Direct unit tests for :class:`GaussianpKaJob` and friends in
``chemsmart.jobs.gaussian.pka``.

CLI-level invocation is covered by ``tests/test_pka.py`` through
``CliRunner``, but that suite mostly exercises backend/batch selection
rather than the job class's own orchestration logic. These tests
exercise the class directly: job wiring (with/without a reference
acid), phase completion tracking, the ``_run`` phase state machine,
output properties, and thermochemistry delegation.

``GaussianOptJob`` and ``GaussianSinglePointJob`` are mocked throughout
so tests don't depend on real Gaussian settings/molecule validation in
those classes -- only on how ``GaussianpKaJob`` wires them together.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.gaussian.pka import (
    GaussianpKaAnalyzeJob,
    GaussianpKaJob,
    GaussianpKaThermoJob,
)
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings


class _FakePKaSettings(GaussianpKaJobSettings):
    """Settings double with deterministic sub-job wiring.

    Overrides the pair/settings-generation methods to return plain
    sentinel strings instead of doing real molecule-copying or
    basis/functional resolution, since ``GaussianOptJob`` and
    ``GaussianSinglePointJob`` are mocked out in these tests and never
    actually inspect these values.
    """

    def conjugate_pair_molecules(self, molecule):
        return "PROT_MOL", "CONJ_MOL"

    def _create_gas_phase_job_settings(self, molecule):
        return "PROT_OPT_SETTINGS", "CONJ_OPT_SETTINGS"

    def reference_pair_molecules(self):
        return "REF_ACID_MOL", "REF_CONJ_MOL"

    def reference_pair_job_settings(self):
        return "REF_ACID_OPT_SETTINGS", "REF_CONJ_OPT_SETTINGS"

    def _create_solution_phase_sp_settings(self, molecule):
        return "PROT_SP_SETTINGS", "CONJ_SP_SETTINGS"

    def reference_pair_sp_job_settings(self):
        return "REF_ACID_SP_SETTINGS", "REF_CONJ_SP_SETTINGS"


@pytest.fixture(autouse=True)
def _clear_shared_reference_cache():
    # `_shared_reference_molecule_cache` is a class attribute shared by
    # every GaussianpKaJob instance; reset it so tests don't leak state.
    GaussianpKaJob._shared_reference_molecule_cache.clear()
    yield
    GaussianpKaJob._shared_reference_molecule_cache.clear()


@pytest.fixture()
def pka_settings():
    return _FakePKaSettings()


@pytest.fixture()
def pka_settings_with_reference():
    return _FakePKaSettings(
        reference_file="href.xyz",
        reference_proton_index=1,
        reference_charge=0,
        reference_multiplicity=1,
    )


def _make_job_factory():
    """Return a side_effect factory for GaussianOptJob/GaussianSinglePointJob.

    Records every created mock job in ``factory.created`` (in creation
    order) so tests can reach in and configure ``is_complete``/
    ``_output``/``outputfile`` per job.
    """

    created = []

    def factory(*, molecule, settings, label, jobrunner, skip_completed):
        job = MagicMock()
        job.molecule = molecule
        job.settings = settings
        job.label = label
        job.jobrunner = jobrunner
        job.skip_completed = skip_completed
        job.is_complete.return_value = False
        job._output.return_value = None
        job.outputfile = f"{label}.log"

        # By default, simulate a successful run: once `.run()` is
        # called the job reports itself complete. Tests that need to
        # simulate a job staying incomplete after running should clear
        # this with `job.run.side_effect = None`.
        def _mark_complete_on_run():
            job.is_complete.return_value = True

        job.run.side_effect = _mark_complete_on_run
        created.append(job)
        return job

    factory.created = created
    return factory


@pytest.fixture()
def opt_factory():
    return _make_job_factory()


@pytest.fixture()
def sp_factory():
    return _make_job_factory()


@pytest.fixture()
def patched_sub_jobs(opt_factory, sp_factory):
    with (
        patch(
            "chemsmart.jobs.gaussian.pka.GaussianOptJob",
            side_effect=opt_factory,
        ),
        patch(
            "chemsmart.jobs.gaussian.pka.GaussianSinglePointJob",
            side_effect=sp_factory,
        ),
    ):
        yield opt_factory, sp_factory


class TestInitValidation:
    def test_rejects_non_pka_settings(self, ethanol_molecule):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        with pytest.raises(ValueError, match="GaussianpKaJobSettings"):
            GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=GaussianJobSettings.default(),
                label="pka_test",
            )


class TestPrepareJobsNoReference:
    def test_creates_only_target_opt_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        opt_factory, _ = patched_sub_jobs
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings,
            label="pka_test",
        )
        assert job.has_reference_jobs is False
        assert job.ref_opt_jobs == []
        assert job.ref_acid_job is None
        assert job.ref_conjugate_base_job is None
        assert len(opt_factory.created) == 2
        assert job.protonated_job.label == "pka_test_HA_opt"
        assert job.conjugate_base_job.label == "pka_test_A_opt"
        assert job.opt_jobs == [job.protonated_job, job.conjugate_base_job]


class TestPrepareJobsWithReference:
    def test_creates_reference_opt_jobs(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        opt_factory, _ = patched_sub_jobs
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.has_reference_jobs is True
        assert len(opt_factory.created) == 4
        assert job.ref_acid_job.label == "pka_test_HRef_opt"
        assert job.ref_conjugate_base_job.label == "pka_test_Ref_opt"
        assert job.ref_opt_jobs == [
            job.ref_acid_job,
            job.ref_conjugate_base_job,
        ]

    def test_reference_pair_molecules_cached_across_jobs(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        with patch.object(
            _FakePKaSettings,
            "reference_pair_molecules",
            return_value=("REF_ACID_MOL", "REF_CONJ_MOL"),
        ) as mocked:
            GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings_with_reference,
                label="pka_job_1",
            )
            GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings_with_reference.copy(),
                label="pka_job_2",
            )
        # Same cache key (scheme, reference_file, indices/charges) ->
        # the expensive geometry-building call only happens once.
        assert mocked.call_count == 1


class TestReferenceCacheHelpers:
    def test_reference_cache_key_none_for_none_settings(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._reference_cache_key(None) is None

    def test_reference_cache_key_none_without_reference_file(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._reference_cache_key(pka_settings) is None

    def test_get_cached_reference_pair_none_without_reference_file(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._get_cached_reference_pair(pka_settings) is None

    def test_prepare_pka_jobs_noop_when_settings_none(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.settings = None
        job._prepare_pka_jobs()  # must not raise

    def test_run_ref_opt_jobs_noop_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_ref_opt_jobs()  # must not raise; nothing to run

    def test_create_ref_sp_jobs_falls_back_when_cache_miss(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        with patch.object(
            job, "_get_cached_reference_pair", return_value=None
        ):
            job._create_ref_sp_jobs()
        assert job.ref_sp_jobs is not None

    def test_prepare_ref_opt_jobs_falls_back_when_cache_miss(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        with patch.object(
            job, "_get_cached_reference_pair", return_value=None
        ):
            ref_acid_job, ref_conjugate_base_job = job._prepare_ref_opt_jobs()
        assert ref_acid_job.label == "pka_test_HRef_opt"
        assert ref_conjugate_base_job.label == "pka_test_Ref_opt"


class TestOptimizedMoleculeFromJob:
    def test_returns_output_molecule_on_normal_termination(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        fake_output = MagicMock(normal_termination=True, molecule="OPT_MOL")
        child_job = MagicMock()
        child_job._output.return_value = fake_output
        result = job._optimized_molecule_from_job(child_job, "FALLBACK")
        assert result == "OPT_MOL"

    def test_falls_back_when_not_normally_terminated(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        fake_output = MagicMock(normal_termination=False, molecule="OPT_MOL")
        child_job = MagicMock()
        child_job._output.return_value = fake_output
        result = job._optimized_molecule_from_job(child_job, "FALLBACK")
        assert result == "FALLBACK"

    def test_falls_back_when_no_output(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        child_job = MagicMock()
        child_job._output.return_value = None
        result = job._optimized_molecule_from_job(child_job, "FALLBACK")
        assert result == "FALLBACK"


class TestCompletionTracking:
    def test_opt_jobs_are_complete_false_when_empty(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.opt_jobs = []
        assert job._opt_jobs_are_complete() is False

    def test_opt_jobs_are_complete_true_when_all_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        assert job._opt_jobs_are_complete() is True

    def test_ref_opt_jobs_are_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._ref_opt_jobs_are_complete() is True

    def test_ref_opt_jobs_are_complete_false_when_incomplete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job._ref_opt_jobs_are_complete() is False
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        assert job._ref_opt_jobs_are_complete() is True

    def test_sp_jobs_are_complete_false_when_none(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._sp_jobs_are_complete() is False

    def test_ref_sp_jobs_are_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._ref_sp_jobs_are_complete() is True

    def test_is_complete_requires_all_phases(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
        sp_factory,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.is_complete() is False

        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        assert job.is_complete() is False  # SP jobs not even created yet

        job._create_sp_jobs()
        for sp_job in job.sp_jobs:
            sp_job.is_complete.return_value = True
        assert job.is_complete() is False  # reference opt still incomplete

        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        assert job.is_complete() is False  # reference SP not created yet

        job._create_ref_sp_jobs()
        for ref_sp_job in job.ref_sp_jobs:
            ref_sp_job.is_complete.return_value = True
        assert job.is_complete() is True

    def test_is_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._create_sp_jobs()
        for sp_job in job.sp_jobs:
            sp_job.is_complete.return_value = True
        assert job.is_complete() is True


class TestRunOptAndSpJobs:
    def test_run_sp_jobs_warns_and_skips_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_sp_jobs()
        assert job.sp_jobs is None

    def test_run_sp_jobs_creates_and_runs_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, sp_factory
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._run_sp_jobs()
        assert job.sp_jobs is not None
        for sp_job in job.sp_jobs:
            sp_job.run.assert_called_once()

    def test_run_ref_sp_jobs_noop_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_ref_sp_jobs()
        assert job.ref_sp_jobs is None

    def test_run_ref_sp_jobs_warns_when_ref_opt_incomplete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job._run_ref_sp_jobs()
        assert job.ref_sp_jobs is None

    def test_run_ref_sp_jobs_creates_and_runs_when_ref_opt_complete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._run_ref_sp_jobs()
        assert job.ref_sp_jobs is not None
        for ref_sp_job in job.ref_sp_jobs:
            ref_sp_job.run.assert_called_once()

    def test_run_sp_jobs_reuses_existing_sp_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._run_sp_jobs()
        existing_sp_jobs = job.sp_jobs

        job._run_sp_jobs()  # second call must not recreate sp_jobs
        assert job.sp_jobs is existing_sp_jobs
        for sp_job in existing_sp_jobs:
            assert sp_job.run.call_count == 2

    def test_run_ref_sp_jobs_reuses_existing_ref_sp_jobs(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._run_ref_sp_jobs()
        existing_ref_sp_jobs = job.ref_sp_jobs

        job._run_ref_sp_jobs()  # second call must not recreate ref_sp_jobs
        assert job.ref_sp_jobs is existing_ref_sp_jobs
        for ref_sp_job in existing_ref_sp_jobs:
            assert ref_sp_job.run.call_count == 2


class TestRunPhaseStateMachine:
    def test_stops_after_incomplete_opt_phase(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        # Both opt jobs stay incomplete after running (clear the
        # factory's default auto-complete-on-run behaviour).
        for opt_job in job.opt_jobs:
            opt_job.run.side_effect = None
        job._run()
        # `_execute_phase_jobs` stops the loop as soon as one job is
        # found incomplete, so only the first job actually ran.
        job.opt_jobs[0].run.assert_called_once()
        job.opt_jobs[1].run.assert_not_called()
        assert job.sp_jobs is None

    def test_runs_sp_phase_when_opt_completes(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        # Default factory behaviour: each job completes once run().
        job._run()
        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        assert job.sp_jobs is not None
        for sp_job in job.sp_jobs:
            sp_job.run.assert_called_once()

    def test_reference_opt_failure_blocks_target_sp_phase(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        # Target opt jobs succeed (default behaviour); reference opt
        # jobs are run but never mark themselves complete.
        for ref_job in job.ref_opt_jobs:
            ref_job.run.side_effect = None
        job._run()
        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        job.ref_opt_jobs[0].run.assert_called_once()
        # Target SP jobs must not have been created/run: reference opt
        # incompleteness halts the whole workflow before that phase.
        assert job.sp_jobs is None

    def test_full_run_with_reference_completes_all_phases(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job._run()
        assert job.is_complete() is True

    def test_stops_when_sp_phase_fails(self, ethanol_molecule, pka_settings):
        opt_factory = _make_job_factory()

        def _sp_factory_never_completes(
            *, molecule, settings, label, jobrunner, skip_completed
        ):
            job = MagicMock()
            job.label = label
            job.outputfile = f"{label}.log"
            job.is_complete.return_value = False
            job._output.return_value = None
            return job

        with (
            patch(
                "chemsmart.jobs.gaussian.pka.GaussianOptJob",
                side_effect=opt_factory,
            ),
            patch(
                "chemsmart.jobs.gaussian.pka.GaussianSinglePointJob",
                side_effect=_sp_factory_never_completes,
            ),
        ):
            job = GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings,
                label="pka_test",
            )
            job._run()

        assert job.sp_jobs is not None
        # `_execute_phase_jobs` stops as soon as one job is found
        # incomplete, so only the first SP job actually ran.
        job.sp_jobs[0].run.assert_called_once()
        job.sp_jobs[1].run.assert_not_called()
        assert job._sp_jobs_are_complete() is False


class TestOutputProperties:
    def test_target_output_properties_delegate_to_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.protonated_job._output.return_value = "HA_OUT"
        job.conjugate_base_job._output.return_value = "A_OUT"
        assert job.protonated_output == "HA_OUT"
        assert job.conjugate_base_output == "A_OUT"

    def test_reference_output_properties_none_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job.ref_acid_output is None
        assert job.ref_conjugate_base_output is None
        assert job.ref_acid_sp_output is None
        assert job.ref_conjugate_base_sp_output is None

    def test_reference_output_properties_delegate_when_present(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job.ref_acid_job._output.return_value = "HREF_OUT"
        job.ref_conjugate_base_job._output.return_value = "REF_OUT"
        assert job.ref_acid_output == "HREF_OUT"
        assert job.ref_conjugate_base_output == "REF_OUT"

    def test_target_sp_output_properties_delegate_to_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._create_sp_jobs()
        job.protonated_sp_job._output.return_value = "HA_SP_OUT"
        job.conjugate_base_sp_job._output.return_value = "A_SP_OUT"
        assert job.protonated_sp_output == "HA_SP_OUT"
        assert job.conjugate_base_sp_output == "A_SP_OUT"

    def test_reference_sp_output_properties_delegate_when_present(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._create_ref_sp_jobs()
        job.ref_acid_sp_job._output.return_value = "HREF_SP_OUT"
        job.ref_conjugate_base_sp_job._output.return_value = "REF_SP_OUT"
        assert job.ref_acid_sp_output == "HREF_SP_OUT"
        assert job.ref_conjugate_base_sp_output == "REF_SP_OUT"


class TestCompletionCheckEdgeCases:
    def test_ref_opt_jobs_are_complete_false_when_list_empty(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job.ref_opt_jobs = []
        assert job._ref_opt_jobs_are_complete() is False

    def test_ref_sp_jobs_are_complete_false_when_none(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.ref_sp_jobs is None
        assert job._ref_sp_jobs_are_complete() is False


class TestPkaOutputFiles:
    def test_excludes_reference_paths_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        ha_file, a_file, href_file, ref_file = job._pka_output_files()
        assert ha_file == "pka_test_HA_opt.log"
        assert a_file == "pka_test_A_opt.log"
        assert href_file is None
        assert ref_file is None

    def test_excludes_reference_paths_when_ref_opt_incomplete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        _, _, href_file, ref_file = job._pka_output_files()
        assert href_file is None
        assert ref_file is None

    def test_includes_reference_paths_when_ref_opt_complete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        _, _, href_file, ref_file = job._pka_output_files()
        assert href_file == "pka_test_HRef_opt.log"
        assert ref_file == "pka_test_Ref_opt.log"


class TestThermochemistryDelegation:
    def test_compute_thermochemistry_raises_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        with pytest.raises(ValueError, match="not complete"):
            job.compute_thermochemistry()

    def test_compute_thermochemistry_delegates_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16pKaOutput"
            ".compute_pka_thermochemistry",
            return_value={"pKa": 4.2},
        ) as mocked:
            result = job.compute_thermochemistry()

        assert result == {"pKa": 4.2}
        mocked.assert_called_once()
        _, kwargs = mocked.call_args
        assert kwargs["ha_file"] == "pka_test_HA_opt.log"
        assert kwargs["a_file"] == "pka_test_A_opt.log"
        assert kwargs["href_file"] is None
        assert kwargs["temperature"] == pka_settings.temperature

    def test_print_thermochemistry_raises_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        with pytest.raises(ValueError, match="not complete"):
            job.print_thermochemistry()

    def test_print_thermochemistry_delegates_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, sp_factory
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._create_sp_jobs()

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16pKaOutput"
            ".print_pka_summary"
        ) as mocked:
            job.print_thermochemistry()

        mocked.assert_called_once()
        _, kwargs = mocked.call_args
        assert kwargs["ha_gas_file"] == "pka_test_HA_opt.log"
        assert kwargs["ha_solv_file"] == "pka_test_HA_sp.log"
        assert kwargs["href_solv_file"] is None

    def test_print_thermochemistry_includes_reference_solv_files(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._create_sp_jobs()
        job._create_ref_sp_jobs()

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16pKaOutput"
            ".print_pka_summary"
        ) as mocked:
            job.print_thermochemistry()

        _, kwargs = mocked.call_args
        assert kwargs["href_solv_file"] == "pka_test_HRef_sp.log"
        assert kwargs["ref_solv_file"] == "pka_test_Ref_sp.log"


class TestMoleculeProperties:
    def test_original_mol_is_input_molecule(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job.original_mol is job.molecule

    def test_conjugate_base_mol_from_job(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.conjugate_base_job.molecule = "CONJ_FROM_JOB"
        assert job.conjugate_base_mol == "CONJ_FROM_JOB"

    def test_conjugate_base_mol_falls_back_without_job(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.conjugate_base_job = None
        assert job.conjugate_base_mol == "CONJ_MOL"


class TestGaussianpKaAnalyzeJob:
    def test_run_calls_print_thermochemistry(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaAnalyzeJob(
            input_file=ethanol_molecule,
            settings=pka_settings,
            label="pka_analyze",
        )
        with patch.object(job, "print_thermochemistry") as mocked:
            job._run()
        mocked.assert_called_once()

    def test_run_swallows_exceptions(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaAnalyzeJob(
            input_file=ethanol_molecule,
            settings=pka_settings,
            label="pka_analyze",
        )
        with patch.object(
            job, "print_thermochemistry", side_effect=RuntimeError("boom")
        ):
            job._run()  # must not raise


class TestGaussianpKaThermoJob:
    def test_type_identifier(self):
        assert GaussianpKaThermoJob.TYPE == "g16pka_thermo"
        assert issubclass(GaussianpKaThermoJob, GaussianpKaAnalyzeJob)
