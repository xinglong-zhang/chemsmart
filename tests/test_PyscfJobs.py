"""Tests for PySCF job objects and the fake job runner.

The fake runner writes the real driver script and a synthetic results file
without importing pyscf, so these tests exercise the
settings -> writer -> output -> parse path in an environment where PySCF is
not installed.
"""

import os

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.hess import PySCFHessJob
from chemsmart.jobs.pyscf.opt import PySCFOptJob
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.pyscf.singlepoint import PySCFSinglePointJob


@pytest.fixture()
def water_molecule(pyscf_inputs_directory):
    filepath = os.path.join(pyscf_inputs_directory, "water.xyz")
    molecule = Molecule.from_filepath(filepath)
    molecule.charge = 0
    molecule.multiplicity = 1
    return molecule


def _settings(jobtype, **extra):
    return PySCFJobSettings(
        functional="b3lyp",
        basis="def2-svp",
        charge=0,
        multiplicity=1,
        jobtype=jobtype,
        **extra,
    )


class TestPySCFJobConstruction:
    """Job identity, file paths, and constructor guards."""

    def test_job_file_paths_follow_the_label(
        self, water_molecule, temporary_working_dir
    ):
        job = PySCFSinglePointJob(
            molecule=water_molecule,
            settings=_settings("sp"),
            label="water_sp",
        )
        folder = str(temporary_working_dir)
        assert job.PROGRAM == "PySCF"
        assert job.inputfile == os.path.join(folder, "water_sp.py")
        assert job.outputfile == os.path.join(folder, "water_sp.out")
        assert job.resultsfile == os.path.join(folder, "water_sp.h5")

    def test_job_without_results_is_not_complete(
        self, water_molecule, temporary_working_dir
    ):
        job = PySCFSinglePointJob(
            molecule=water_molecule,
            settings=_settings("sp"),
            label="water_sp",
        )
        assert not os.path.exists(job.resultsfile)
        assert job.is_complete() is False

    def test_settings_must_be_pyscf_settings(
        self, water_molecule, temporary_working_dir
    ):
        with pytest.raises(ValueError, match="Settings must be instance of"):
            PySCFSinglePointJob(
                molecule=water_molecule,
                settings={"functional": "b3lyp"},
                label="water_sp",
            )

    def test_molecule_must_be_a_molecule(self, temporary_working_dir):
        with pytest.raises(ValueError):
            PySCFSinglePointJob(
                molecule="water",
                settings=_settings("sp"),
                label="water_sp",
            )


class TestFakePySCFJobRunner:
    """The preview runner must produce a parseable, clearly fake result."""

    @staticmethod
    def _run(job_class, jobtype, molecule, jobrunner, **extra):
        job = job_class(
            molecule=molecule,
            settings=_settings(jobtype, **extra),
            label=f"water_{jobtype}",
            jobrunner=jobrunner,
        )
        assert jobrunner.run(job) == 0
        return job

    def test_preview_writes_script_log_and_results(
        self,
        water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            water_molecule,
            pyscf_jobrunner_no_scratch,
        )
        assert os.path.exists(job.inputfile)
        assert os.path.exists(job.outputfile)
        assert os.path.exists(job.resultsfile)

    def test_preview_result_is_marked_fake(
        self,
        water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            water_molecule,
            pyscf_jobrunner_no_scratch,
        )
        output = PySCFOutput(job.outputfile)
        # A preview must never be mistaken for a calculation: the committed
        # real artifacts report the producing PySCF version instead.
        assert output.version == "fake"
        assert output.engine == "cpu"

    @pytest.mark.parametrize(
        ("job_class", "jobtype", "extra"),
        [
            (PySCFOptJob, "opt", {}),
            (PySCFHessJob, "hess", {"freq": True}),
        ],
    )
    def test_preview_generates_a_driver_script_per_jobtype(
        self,
        job_class,
        jobtype,
        extra,
        water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            job_class,
            jobtype,
            water_molecule,
            pyscf_jobrunner_no_scratch,
            **extra,
        )
        with open(job.inputfile) as handle:
            script = handle.read()
        assert "ChemSmart PySCF driver" in script
        assert "def2-svp" in script

    def test_preview_in_scratch_copies_results_back(
        self,
        water_molecule,
        pyscf_jobrunner_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            water_molecule,
            pyscf_jobrunner_scratch,
        )
        # The calculation runs in the scratch directory, but the artifacts
        # must end up next to the job.
        assert os.path.exists(job.resultsfile)
        assert os.path.exists(job.outputfile)
        assert PySCFOutput(job.outputfile).version == "fake"

    def test_preview_records_the_requested_solvent(
        self,
        water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            water_molecule,
            pyscf_jobrunner_no_scratch,
            solvent_model="cpcm",
            solvent_id="water",
        )
        output = PySCFOutput(job.outputfile)
        assert output.solvent_on
        assert output.solvent_model == "cpcm"
        assert output.solvent_id == "water"
