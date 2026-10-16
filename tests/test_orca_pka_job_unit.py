"""
Direct unit tests for :class:`ORCApKaJob` in ``chemsmart.jobs.orca.pka``.

CLI-level invocation is covered by ``tests/test_pka.py`` through
``CliRunner``, but that suite mostly exercises backend/batch selection
rather than the job class's own orchestration logic. These tests
exercise the class directly: proton_index validation, lazily-cached
job-list properties, the file-backed completion check
(``_bind_subjob``/``_subjob_is_complete``/``_subjob_output``), the
``_run`` phase state machine, output properties, and thermochemistry
delegation.

``ORCAOptJob`` and ``ORCASinglePointJob`` are mocked throughout so
tests don't depend on real ORCA settings/molecule validation in those
classes -- only on how ``ORCApKaJob`` wires them together. Real
temporary ``.out`` files are used only for the dedicated
``_subjob_is_complete``/``_subjob_output`` tests, since those methods
read from disk.
"""

from unittest.mock import MagicMock, PropertyMock, patch

import pytest

from chemsmart.jobs.orca.pka import ORCApKaJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings


class _FakePKaSettings(ORCApKaJobSettings):
    """Settings double with deterministic sub-job wiring.

    Overrides the pair/settings-generation methods to return plain
    sentinel strings instead of doing real molecule-copying or
    basis/functional resolution, since ``ORCAOptJob`` and
    ``ORCASinglePointJob`` are mocked out in these tests and never
    actually inspect these values.
    """

    def conjugate_pair_molecules(self, molecule):
        return "PROT_MOL", "CONJ_MOL"

    def conjugate_pair_job_settings(self, molecule):
        return "PROT_OPT_SETTINGS", "CONJ_OPT_SETTINGS"

    def reference_pair_molecules(self):
        return "REF_ACID_MOL", "REF_CONJ_MOL"

    def reference_pair_job_settings(self):
        return "REF_ACID_OPT_SETTINGS", "REF_CONJ_OPT_SETTINGS"

    def _create_solution_phase_sp_settings(self, molecule):
        return "PROT_SP_SETTINGS", "CONJ_SP_SETTINGS"

    def reference_pair_sp_job_settings(self):
        return "REF_ACID_SP_SETTINGS", "REF_CONJ_SP_SETTINGS"

    def get_reference_molecule(self):
        return "REF_ACID_MOL_FALLBACK"

    def get_reference_conjugate_base_molecule(self):
        return "REF_CONJ_MOL_FALLBACK"


@pytest.fixture(autouse=True)
def _clear_shared_reference_cache():
    ORCApKaJob._shared_reference_molecule_cache.clear()
    yield
    ORCApKaJob._shared_reference_molecule_cache.clear()


@pytest.fixture()
def pka_settings():
    return _FakePKaSettings(proton_index=1)


@pytest.fixture()
def pka_settings_with_reference():
    return _FakePKaSettings(
        proton_index=1,
        reference_file="href.xyz",
        reference_proton_index=1,
        reference_charge=0,
        reference_multiplicity=1,
    )


def _make_job_factory():
    """Return a side_effect factory for ORCAOptJob/ORCASinglePointJob.

    Records every created mock job (in creation order) so tests can
    reach in and configure them. Note that ``_bind_subjob`` always
    overwrites ``job.is_complete`` right after creation in production
    code, so orchestration tests must (re)configure ``is_complete``
    *after* pulling jobs out via the job-list properties -- see
    ``_stub_is_complete`` below.
    """

    created = []

    def factory(*, molecule, settings, label, jobrunner, skip_completed):
        job = MagicMock()
        job.molecule = molecule
        job.settings = settings
        job.label = label
        job.jobrunner = jobrunner
        job.skip_completed = skip_completed
        job._output.return_value = None
        job.outputfile = f"{label}.out"
        created.append(job)
        return job

    factory.created = created
    return factory


def _stub_is_complete(job, auto_complete_on_run=True):
    """(Re-)configure a mock job's completion behaviour.

    Must be called *after* the job has been produced by one of the
    lazy job-list properties, since ``_bind_subjob`` unconditionally
    replaces ``job.is_complete`` with a real file-checking closure as
    part of job creation.
    """
    job.is_complete = MagicMock(return_value=False)
    if auto_complete_on_run:

        def _mark_complete():
            job.is_complete.return_value = True

        job.run.side_effect = _mark_complete
    else:
        job.run.side_effect = None
    return job


def _disable_bind_subjob(job):
    """Replace real file-backed completion binding with a plain mock.

    ``_run_sp_jobs``/``_run_ref_sp_jobs`` reset ``_sp_jobs``/
    ``_ref_sp_jobs`` to ``None`` before creating fresh SP jobs (so they
    pick up the latest optimized geometry), which means any earlier
    stubbing done on a previously-fetched ``sp_jobs`` list is
    discarded, and the freshly created jobs get real
    ``_subjob_is_complete`` bindings pointed at output files that don't
    exist on disk. Call this right after constructing the job (before
    touching any job-list property) to make every subsequently-created
    sub-job's ``is_complete`` a plain, freely-configurable mock
    instead -- the real binding is covered separately by
    ``TestSubjobOutputHelpers``.
    """

    def _stub_bind(sub_job, legacy_label=None):
        _stub_is_complete(sub_job)

    job._bind_subjob = _stub_bind


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
            "chemsmart.jobs.orca.pka.ORCAOptJob",
            side_effect=opt_factory,
        ),
        patch(
            "chemsmart.jobs.orca.pka.ORCASinglePointJob",
            side_effect=sp_factory,
        ),
    ):
        yield opt_factory, sp_factory


class TestInitValidation:
    def test_rejects_non_pka_settings(self, ethanol_molecule):
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        with pytest.raises(ValueError, match="ORCApKaJobSettings"):
            ORCApKaJob(
                molecule=ethanol_molecule,
                settings=ORCAJobSettings.default(),
                label="pka_test",
            )

    def test_rejects_missing_proton_index(self, ethanol_molecule):
        with pytest.raises(ValueError, match="proton_index"):
            ORCApKaJob(
                molecule=ethanol_molecule,
                settings=_FakePKaSettings(proton_index=None),
                label="pka_test",
            )


class TestLazyJobProperties:
    def test_opt_jobs_created_lazily_and_cached(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        opt_factory, _ = patched_sub_jobs
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert opt_factory.created == []  # not created until accessed
        opt_jobs = job.opt_jobs
        assert len(opt_factory.created) == 2
        assert job.opt_jobs is opt_jobs  # cached, not recreated
        assert job.protonated_job.label == "pka_test_HA_opt"
        assert job.conjugate_base_job.label == "pka_test_A_opt"

    def test_sp_jobs_created_lazily_and_cached(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        _, sp_factory = patched_sub_jobs
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        sp_jobs = job.sp_jobs
        assert len(sp_factory.created) == 2
        assert job.sp_jobs is sp_jobs
        assert job.protonated_sp_job.label == "pka_test_HA_sp"
        assert job.conjugate_base_sp_job.label == "pka_test_A_sp"

    def test_ref_job_properties_none_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job.has_reference_jobs is False
        assert job.ref_opt_jobs is None
        assert job.ref_acid_job is None
        assert job.ref_conjugate_base_job is None
        assert job.ref_sp_jobs is None
        assert job.ref_acid_sp_job is None
        assert job.ref_conjugate_base_sp_job is None

    def test_ref_opt_jobs_created_lazily_with_reference(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.has_reference_jobs is True
        ref_opt_jobs = job.ref_opt_jobs
        assert job.ref_opt_jobs is ref_opt_jobs
        assert job.ref_acid_job.label == "href"
        assert job.ref_conjugate_base_job.label == "href_cb"

    def test_ref_sp_jobs_created_lazily_with_reference(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        ref_sp_jobs = job.ref_sp_jobs
        assert job.ref_sp_jobs is ref_sp_jobs
        assert job.ref_acid_sp_job.label == "href_sp"
        assert job.ref_conjugate_base_sp_job.label == "href_cb_sp"

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
            job1 = ORCApKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings_with_reference,
                label="pka_job_1",
            )
            _ = job1.ref_opt_jobs
            job2 = ORCApKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings_with_reference.copy(),
                label="pka_job_2",
            )
            _ = job2.ref_opt_jobs
        assert mocked.call_count == 1

    def test_reference_molecule_properties_use_cache(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.reference_molecule == "REF_ACID_MOL"
        assert job.reference_conjugate_base_molecule == "REF_CONJ_MOL"

    def test_reference_molecule_properties_fall_back_without_cache(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        with patch.object(
            job, "_get_cached_reference_pair", return_value=None
        ):
            assert job.reference_molecule == "REF_ACID_MOL_FALLBACK"
            assert (
                job.reference_conjugate_base_molecule
                == "REF_CONJ_MOL_FALLBACK"
            )

    def test_prepare_ref_opt_jobs_falls_back_when_cache_miss(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        with patch.object(
            job, "_get_cached_reference_pair", return_value=None
        ):
            ref_acid_job, ref_cb_job = job._prepare_ref_opt_jobs()
        assert ref_acid_job.label == "href"
        assert ref_cb_job.label == "href_cb"


class TestReferenceCacheHelpers:
    def test_reference_cache_key_none_for_none_settings(self):
        assert ORCApKaJob._reference_cache_key(None) is None

    def test_reference_cache_key_none_without_reference_file(
        self, pka_settings
    ):
        assert ORCApKaJob._reference_cache_key(pka_settings) is None

    def test_get_cached_reference_pair_none_without_reference_file(
        self, pka_settings
    ):
        assert ORCApKaJob._get_cached_reference_pair(pka_settings) is None


class TestRefBasenameHelpers:
    def test_ref_basename_none_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._ref_basename is None
        assert job._ref_conjugate_base_label is None

    def test_ref_basename_derived_from_reference_file(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job._ref_basename == "href"
        assert job._ref_conjugate_base_label == "href_cb"


class TestMoleculeProperties:
    def test_protonated_and_conjugate_base_molecule(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job.protonated_molecule == "PROT_MOL"
        assert job.conjugate_base_molecule == "CONJ_MOL"

    def test_settings_class(self):
        assert ORCApKaJob.settings_class() is ORCApKaJobSettings

    def test_prepare_sp_jobs_uses_optimized_geometry_when_available(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        protonated_job = job.protonated_job
        conjugate_base_job = job.conjugate_base_job

        def fake_subjob_output(sub_job, legacy_label=None):
            if sub_job is protonated_job:
                return MagicMock(molecule="HA_OPT_MOL")
            if sub_job is conjugate_base_job:
                return MagicMock(molecule="A_OPT_MOL")
            return None

        with patch.object(
            job, "_subjob_output", side_effect=fake_subjob_output
        ):
            sp_jobs = job._prepare_sp_jobs()

        assert sp_jobs[0].molecule == "HA_OPT_MOL"
        assert sp_jobs[1].molecule == "A_OPT_MOL"

    def test_prepare_ref_sp_jobs_uses_optimized_geometry_when_available(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        ref_acid_job = job.ref_acid_job
        ref_cb_job = job.ref_conjugate_base_job
        ref_acid_job._output.return_value = MagicMock(
            normal_termination=True, molecule="HREF_OPT_MOL"
        )
        ref_cb_job._output.return_value = MagicMock(
            normal_termination=True, molecule="REF_OPT_MOL"
        )

        ref_sp_jobs = job._prepare_ref_sp_jobs()
        assert ref_sp_jobs[0].molecule == "HREF_OPT_MOL"
        assert ref_sp_jobs[1].molecule == "REF_OPT_MOL"

    def test_prepare_ref_jobs_bind_without_legacy_label_when_no_basename(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        with patch.object(
            type(job),
            "_ref_basename",
            new_callable=PropertyMock,
            return_value=None,
        ):
            ref_acid_job, ref_cb_job = job._prepare_ref_opt_jobs()
            ref_acid_sp_job, ref_cb_sp_job = job._prepare_ref_sp_jobs()
        assert ref_acid_job is not None
        assert ref_acid_sp_job is not None


class TestMakeSpJob:
    """`_make_sp_job` is defined but not called anywhere in production
    code (dead helper); tested directly for both geometry branches."""

    def test_uses_optimized_geometry_on_normal_termination(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        opt_job = MagicMock()
        opt_job._output.return_value = MagicMock(
            normal_termination=True, molecule="OPT_MOL"
        )
        sp_job = job._make_sp_job(
            opt_job, "FALLBACK_MOL", "SP_SETTINGS", "sp_label"
        )
        assert sp_job.molecule == "OPT_MOL"

    def test_falls_back_when_not_normally_terminated(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        opt_job = MagicMock()
        opt_job._output.return_value = MagicMock(
            normal_termination=False, molecule="OPT_MOL"
        )
        sp_job = job._make_sp_job(
            opt_job, "FALLBACK_MOL", "SP_SETTINGS", "sp_label"
        )
        assert sp_job.molecule == "FALLBACK_MOL"


class TestSubjobOutputHelpers:
    """Tests for the file-backed completion-checking helpers.

    Unlike the rest of this file, these use real temporary ``.out``
    files and patch ``ORCAOutput`` (imported locally inside these
    methods) rather than mocking ORCAOptJob/ORCASinglePointJob.
    """

    def _make_job(self, tmp_path, folder, outputfile):
        job = MagicMock()
        job.jobrunner = None
        job.outputfile = str(tmp_path / outputfile)
        return job

    def test_output_paths_include_runner_output_and_legacy(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = "job.out"
        child.jobrunner = MagicMock(job_outputfile="scratch/job.out")
        paths = job._subjob_output_paths(child, legacy_label="legacy")
        assert paths == [
            "scratch/job.out",
            "job.out",
            str(tmp_path / "legacy.out"),
        ]

    def test_output_paths_dedupes_and_skips_none_runner(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = "job.out"
        child.jobrunner = None
        paths = job._subjob_output_paths(child, legacy_label=None)
        assert paths == ["job.out"]

    def test_output_paths_skips_falsy_runner_outputfile(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = "job.out"
        child.jobrunner = MagicMock(job_outputfile=None)
        paths = job._subjob_output_paths(child, legacy_label=None)
        assert paths == ["job.out"]

    def test_output_paths_dedupes_duplicate_candidates(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = "job.out"
        child.jobrunner = MagicMock(job_outputfile="job.out")
        paths = job._subjob_output_paths(child, legacy_label=None)
        assert paths == ["job.out"]

    def test_subjob_is_complete_continues_to_next_candidate(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        first = tmp_path / "first.out"
        first.write_text("x")
        second = tmp_path / "legacy.out"
        second.write_text("x")
        child = MagicMock()
        child.outputfile = str(first)
        child.jobrunner = None
        outputs = {
            str(first): MagicMock(normal_termination=False),
            str(second): MagicMock(normal_termination=True),
        }
        with patch(
            "chemsmart.io.orca.output.ORCAOutput",
            side_effect=lambda path: outputs[path],
        ):
            assert (
                job._subjob_is_complete(child, legacy_label="legacy") is True
            )

    def test_subjob_output_continues_to_next_candidate(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        first = tmp_path / "first.out"
        first.write_text("x")
        second = tmp_path / "legacy.out"
        second.write_text("x")
        child = MagicMock()
        child.outputfile = str(first)
        child.jobrunner = None
        second_output = MagicMock(normal_termination=True)
        outputs = {
            str(first): MagicMock(normal_termination=False),
            str(second): second_output,
        }
        with patch(
            "chemsmart.io.orca.output.ORCAOutput",
            side_effect=lambda path: outputs[path],
        ):
            assert (
                job._subjob_output(child, legacy_label="legacy")
                is second_output
            )

    def test_subjob_is_complete_false_when_no_file_exists(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = str(tmp_path / "missing.out")
        child.jobrunner = None
        assert job._subjob_is_complete(child) is False

    def test_subjob_is_complete_true_on_normal_termination(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        outfile = tmp_path / "done.out"
        outfile.write_text("dummy")
        child = MagicMock()
        child.outputfile = str(outfile)
        child.jobrunner = None
        with patch("chemsmart.io.orca.output.ORCAOutput") as mocked_cls:
            mocked_cls.return_value.normal_termination = True
            assert job._subjob_is_complete(child) is True

    def test_subjob_is_complete_false_when_parsing_raises(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        outfile = tmp_path / "broken.out"
        outfile.write_text("dummy")
        child = MagicMock()
        child.outputfile = str(outfile)
        child.jobrunner = None
        with patch(
            "chemsmart.io.orca.output.ORCAOutput",
            side_effect=RuntimeError("bad parse"),
        ):
            assert job._subjob_is_complete(child) is False

    def test_subjob_output_returns_none_when_file_missing(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = str(tmp_path / "missing.out")
        child.jobrunner = None
        assert job._subjob_output(child) is None

    def test_subjob_output_returns_output_on_normal_termination(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        outfile = tmp_path / "done.out"
        outfile.write_text("dummy")
        child = MagicMock()
        child.outputfile = str(outfile)
        child.jobrunner = None
        fake_output = MagicMock(normal_termination=True)
        with patch(
            "chemsmart.io.orca.output.ORCAOutput",
            return_value=fake_output,
        ):
            assert job._subjob_output(child) is fake_output

    def test_subjob_output_skips_files_with_parse_errors(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        outfile = tmp_path / "broken.out"
        outfile.write_text("dummy")
        child = MagicMock()
        child.outputfile = str(outfile)
        child.jobrunner = None
        with patch(
            "chemsmart.io.orca.output.ORCAOutput",
            side_effect=RuntimeError("bad parse"),
        ):
            assert job._subjob_output(child) is None

    def test_bind_subjob_sets_folder_and_is_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, tmp_path
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.folder = str(tmp_path)
        child = MagicMock()
        child.outputfile = str(tmp_path / "missing.out")
        child.jobrunner = None
        job._bind_subjob(child)
        assert child.folder == str(tmp_path)
        assert child.is_complete() is False


class TestRunOrchestration:
    def test_run_stops_after_incomplete_opt_phase(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job, auto_complete_on_run=False)
        job._run()
        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        assert job._sp_jobs is None

    def test_run_completes_all_phases_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        # `_run_sp_jobs` resets `_sp_jobs` to None and recreates fresh
        # SP jobs every call, so pre-stubbing `job.sp_jobs` here
        # wouldn't survive into `_run()` -- disable the real
        # file-backed binding instead (covered separately by
        # TestSubjobOutputHelpers) so every sub-job this job creates
        # auto-completes once run.
        _disable_bind_subjob(job)
        job._run()
        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        assert job.sp_jobs is not None
        for sp_job in job.sp_jobs:
            sp_job.run.assert_called_once()

    def test_run_with_reference_runs_all_four_phases(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        _disable_bind_subjob(job)

        job._run()

        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        for ref_job in job.ref_opt_jobs:
            ref_job.run.assert_called_once()
        assert job.sp_jobs is not None
        for sp_job in job.sp_jobs:
            sp_job.run.assert_called_once()
        assert job.ref_sp_jobs is not None
        for ref_sp_job in job.ref_sp_jobs:
            ref_sp_job.run.assert_called_once()

    def test_run_reference_opt_failure_blocks_sp_phase(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
        for ref_job in job.ref_opt_jobs:
            _stub_is_complete(ref_job, auto_complete_on_run=False)

        job._run()

        for ref_job in job.ref_opt_jobs:
            ref_job.run.assert_called_once()
        assert job._sp_jobs is None

    def test_run_stops_when_sp_phase_fails(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        # `_prepare_opt_jobs` binds the two opt jobs first, then
        # `_prepare_sp_jobs` binds the two SP jobs: make the first two
        # binds auto-complete (opt succeeds) and the rest never
        # complete (SP fails), all through one shared `_bind_subjob`.
        bind_calls = {"n": 0}

        def _stub_bind(sub_job, legacy_label=None):
            bind_calls["n"] += 1
            sub_job.is_complete = MagicMock(return_value=False)
            if bind_calls["n"] <= 2:
                _stub_is_complete(sub_job)  # opt jobs: auto-complete
            else:
                sub_job.run.side_effect = None  # sp jobs: never complete

        job._bind_subjob = _stub_bind

        job._run()

        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        assert job.sp_jobs is not None
        # `stop_on_incomplete=True` breaks the SP loop after the first
        # job is found incomplete.
        job.sp_jobs[0].run.assert_called_once()
        job.sp_jobs[1].run.assert_not_called()

    def test_run_ref_opt_jobs_noop_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_ref_opt_jobs()  # must not raise; nothing to run

    def test_run_ref_sp_jobs_noop_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_ref_sp_jobs()  # must not raise; nothing to run

    def test_run_returns_early_if_sp_jobs_unexpectedly_none(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True

        # `_prepare_sp_jobs` always returns a real tuple in practice,
        # so `self.sp_jobs is None` in `_run` is a defensive branch
        # that can't happen through normal use; force it directly.
        with patch.object(job, "_prepare_sp_jobs", return_value=None):
            job._run()  # must not raise
        assert job._sp_jobs is None


class TestCompletionTracking:
    def test_is_complete_false_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job, auto_complete_on_run=False)
        assert job.is_complete() is False

    def test_is_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True
        for sp_job in job.sp_jobs:
            _stub_is_complete(sp_job)
            sp_job.is_complete.return_value = True
        assert job.is_complete() is True

    def test_is_complete_requires_reference_phases(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True
        for sp_job in job.sp_jobs:
            _stub_is_complete(sp_job)
            sp_job.is_complete.return_value = True
        assert job.is_complete() is False  # ref opt/sp not done

        for ref_job in job.ref_opt_jobs:
            _stub_is_complete(ref_job)
            ref_job.is_complete.return_value = True
        assert job.is_complete() is False  # ref sp not created/done

        for ref_sp_job in job.ref_sp_jobs:
            _stub_is_complete(ref_sp_job)
            ref_sp_job.is_complete.return_value = True
        assert job.is_complete() is True

    def test_opt_jobs_are_complete_helper(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job, auto_complete_on_run=False)
        assert job._opt_jobs_are_complete() is False
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        assert job._opt_jobs_are_complete() is True

    def test_ref_opt_jobs_are_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._ref_opt_jobs_are_complete() is True

    def test_is_complete_false_when_sp_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True
        for sp_job in job.sp_jobs:
            _stub_is_complete(sp_job)  # left incomplete (False)
        assert job.is_complete() is False

    def test_opt_jobs_are_complete_false_when_opt_jobs_empty(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        with patch.object(
            type(job), "opt_jobs", new_callable=PropertyMock, return_value=[]
        ):
            assert job._opt_jobs_are_complete() is False

    def test_ref_opt_jobs_are_complete_false_when_ref_opt_jobs_empty(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        with patch.object(
            type(job),
            "ref_opt_jobs",
            new_callable=PropertyMock,
            return_value=[],
        ):
            assert job._ref_opt_jobs_are_complete() is False


class TestOutputProperties:
    def test_target_output_properties_delegate_to_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.protonated_job._output.return_value = "HA_OUT"
        job.conjugate_base_job._output.return_value = "A_OUT"
        assert job.protonated_output == "HA_OUT"
        assert job.conjugate_base_output == "A_OUT"

    def test_target_sp_output_properties_delegate_to_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.protonated_sp_job._output.return_value = "HA_SP_OUT"
        job.conjugate_base_sp_job._output.return_value = "A_SP_OUT"
        assert job.protonated_sp_output == "HA_SP_OUT"
        assert job.conjugate_base_sp_output == "A_SP_OUT"

    def test_reference_output_properties_none_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
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
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        # Access the lazy ref-SP properties first: creating them reads
        # `ref_acid_job._output()` / `ref_conjugate_base_job._output()`
        # for a possible optimized geometry, so those must still be
        # harmless default mocks (not the sentinel strings below) at
        # creation time.
        ref_acid_sp_job = job.ref_acid_sp_job
        ref_conjugate_base_sp_job = job.ref_conjugate_base_sp_job

        job.ref_acid_job._output.return_value = "HREF_OUT"
        job.ref_conjugate_base_job._output.return_value = "REF_OUT"
        ref_acid_sp_job._output.return_value = "HREF_SP_OUT"
        ref_conjugate_base_sp_job._output.return_value = "REF_SP_OUT"
        assert job.ref_acid_output == "HREF_OUT"
        assert job.ref_conjugate_base_output == "REF_OUT"
        assert job.ref_acid_sp_output == "HREF_SP_OUT"
        assert job.ref_conjugate_base_sp_output == "REF_SP_OUT"


class TestPkaOutputFiles:
    def test_excludes_reference_paths_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        ha_file, a_file, href_file, ref_file = job._pka_output_files()
        assert ha_file == "pka_test_HA_opt.out"
        assert a_file == "pka_test_A_opt.out"
        assert href_file is None
        assert ref_file is None

    def test_includes_reference_paths_when_ref_opt_complete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            _stub_is_complete(ref_job)
            ref_job.is_complete.return_value = True
        _, _, href_file, ref_file = job._pka_output_files()
        assert href_file == "href.out"
        assert ref_file == "href_cb.out"


class TestThermochemistryDelegation:
    def test_compute_thermochemistry_raises_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job, auto_complete_on_run=False)
        with pytest.raises(ValueError, match="not complete"):
            job.compute_thermochemistry()

    def test_compute_thermochemistry_delegates_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True

        with patch(
            "chemsmart.io.orca.output.ORCApKaOutput.compute_pka_thermochemistry",
            return_value={"pKa": 4.2},
        ) as mocked:
            result = job.compute_thermochemistry()

        assert result == {"pKa": 4.2}
        _, kwargs = mocked.call_args
        assert kwargs["ha_file"] == "pka_test_HA_opt.out"
        assert kwargs["href_file"] is None
        assert kwargs["temperature"] == pka_settings.temperature

    def test_print_thermochemistry_raises_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job, auto_complete_on_run=False)
        with pytest.raises(ValueError, match="not complete"):
            job.print_thermochemistry()

    def test_print_thermochemistry_delegates_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True

        with patch(
            "chemsmart.io.orca.output.ORCApKaOutput.print_pka_summary"
        ) as mocked:
            job.print_thermochemistry()

        _, kwargs = mocked.call_args
        assert kwargs["ha_gas_file"] == "pka_test_HA_opt.out"
        assert kwargs["ha_solv_file"] == "pka_test_HA_sp.out"
        assert kwargs["href_solv_file"] is None

    def test_print_thermochemistry_includes_reference_solv_files(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = ORCApKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for opt_job in job.opt_jobs:
            _stub_is_complete(opt_job)
            opt_job.is_complete.return_value = True
        for ref_job in job.ref_opt_jobs:
            _stub_is_complete(ref_job)
            ref_job.is_complete.return_value = True

        with patch(
            "chemsmart.io.orca.output.ORCApKaOutput.print_pka_summary"
        ) as mocked:
            job.print_thermochemistry()

        _, kwargs = mocked.call_args
        assert kwargs["href_solv_file"] == "href_sp.out"
        assert kwargs["ref_solv_file"] == "href_cb_sp.out"
