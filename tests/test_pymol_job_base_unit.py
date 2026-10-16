"""
Direct unit tests for chemsmart.jobs.mol.job.PyMOLJob, the base class
for all PyMOL visualization jobs. No existing test file constructs
PyMOLJob directly -- coverage only came incidentally through its
concrete subclasses -- leaving several of its own branches (explicit
non-None style-default kwargs, _backup_files, _output/_job_is_complete,
and the from_filename/from_pubchem explicit-label/explicit-jobrunner
paths) untested.
"""

import os

import pytest

from chemsmart.jobs.mol.job import PyMOLJob
from chemsmart.jobs.mol.visualize import PyMOLVisualizationJob
from chemsmart.jobs.runner import JobRunner


class TestPyMOLJobInitDefaults:
    def test_defaults_applied_when_values_omitted(self, ethanol_molecule):
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        assert job.isosurface_value == 0.05
        assert job.transparency_value == 0.2
        assert job.surface_quality == 3
        assert job.antialias_value == 3
        assert job.ray_trace_mode == 1

    def test_explicit_values_are_not_overridden(self, ethanol_molecule):
        job = PyMOLJob(
            molecule=ethanol_molecule,
            label="test",
            isosurface_value=0.1,
            transparency_value=0.5,
            surface_quality=1,
            antialias_value=2,
            ray_trace_mode=0,
        )
        assert job.isosurface_value == 0.1
        assert job.transparency_value == 0.5
        assert job.surface_quality == 1
        assert job.antialias_value == 2
        assert job.ray_trace_mode == 0

    def test_source_basename_defaults_to_label(self, ethanol_molecule):
        job = PyMOLJob(molecule=ethanol_molecule, label="mylabel")
        assert job.source_basename == "mylabel"


class TestPyMOLJobFilePaths:
    def test_file_paths_derived_from_label_and_folder(
        self, ethanol_molecule, tmp_path
    ):
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        job.folder = str(tmp_path)
        assert job.inputfile == str(tmp_path / "test.xyz")
        assert job.logfile == str(tmp_path / "log.test")
        assert job.outputfile == str(tmp_path / "test.pse")
        assert job.errfile == str(tmp_path / "test.err")


class TestPyMOLJobOutputAndCompletion:
    def test_output_is_none_when_file_missing(
        self, ethanol_molecule, tmp_path
    ):
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        job.folder = str(tmp_path)
        assert job._output() is None
        assert job._job_is_complete() is False

    def test_output_returns_path_when_file_exists(
        self, ethanol_molecule, tmp_path
    ):
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        job.folder = str(tmp_path)
        (tmp_path / "test.pse").write_text("dummy session")
        assert job._output() == os.path.abspath(job.outputfile)
        assert job._job_is_complete() is True


class TestPyMOLJobBackupFiles:
    def test_backup_files_noop_when_nothing_exists(
        self, ethanol_molecule, tmp_path
    ):
        """backup_file() no-ops for files that don't exist, so this
        just exercises the folder-name/backup_file call sequence."""
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        job.folder = str(tmp_path)
        job._backup_files()  # should not raise

    def test_backup_files_copies_existing_input_and_output(
        self, ethanol_molecule, tmp_path
    ):
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        job.folder = str(tmp_path)
        (tmp_path / "test.xyz").write_text("3\n\nO 0 0 0\n")
        (tmp_path / "test.pse").write_text("dummy session")

        job._backup_files()

        backup_dirs = [
            p
            for p in tmp_path.iterdir()
            if p.is_dir() and p.name.startswith("bk.")
        ]
        assert len(backup_dirs) == 1
        assert (backup_dirs[0] / "test.xyz").exists()
        assert (backup_dirs[0] / "test.pse").exists()

    def test_backup_chk_crashes_since_pymol_jobs_have_no_chkfile(
        self, ethanol_molecule, tmp_path
    ):
        """PyMOLJob (and its subclasses) never define a chkfile
        property -- backup_chk=True unconditionally accesses
        self.chkfile, which doesn't exist for any PyMOL job."""
        job = PyMOLJob(molecule=ethanol_molecule, label="test")
        job.folder = str(tmp_path)
        with pytest.raises(AttributeError, match="chkfile"):
            job._backup_files(backup_chk=True)


class TestPyMOLJobFromFilename:
    """These classmethods are inherited from PyMOLJob but exercised
    here via a concrete subclass, since the base class's own TYPE has
    no matching registered runner for JobRunner.from_job to find."""

    def test_from_filename_uses_explicit_label(
        self, tmp_path, single_molecule_xyz_file
    ):
        job = PyMOLVisualizationJob.from_filename(
            filename=single_molecule_xyz_file,
            label="explicit_label",
            jobrunner=object(),
        )
        assert job.label == "explicit_label"

    def test_from_filename_creates_jobrunner_when_omitted(
        self, single_molecule_xyz_file
    ):
        job = PyMOLVisualizationJob.from_filename(
            filename=single_molecule_xyz_file
        )
        assert isinstance(job.jobrunner, JobRunner)


class TestPyMOLJobFromPubchem:
    def test_from_pubchem_creates_jobrunner_when_omitted(self, monkeypatch):
        from chemsmart.io.molecules.structure import Molecule

        monkeypatch.setattr(
            Molecule,
            "from_pubchem",
            classmethod(lambda cls, identifier: cls_molecule_stub()),
        )
        job = PyMOLVisualizationJob.from_pubchem(
            identifier="water", label="pubchem_test"
        )
        assert isinstance(job.jobrunner, JobRunner)
        assert job.label == "pubchem_test"


def cls_molecule_stub():
    from chemsmart.io.molecules.structure import Molecule

    return Molecule(
        symbols=["O", "H", "H"],
        positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        charge=0,
        multiplicity=1,
    )
