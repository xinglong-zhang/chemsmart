from pathlib import Path

from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.settings import CRESTJobSettings


class TestCRESTJobs:
    def test_conformers_job_with_constraints(
        self, temporary_working_dir, water_molecule, crest_jobrunner_no_scratch
    ):
        settings = CRESTJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="conformers",
            solvent_model="alpb",
            solvent_id="water",
            constraints=[[1, 2]],
            force_constant=0.5,
            gfn_version="gfn2",
        )
        job = CRESTConformerSearchJob(
            molecule=water_molecule,
            settings=settings,
            label="water_conformers",
            jobrunner=crest_jobrunner_no_scratch,
        )
        assert isinstance(job, CRESTConformerSearchJob)
        assert crest_jobrunner_no_scratch.run(job) == 0

        folder = Path(job.folder)
        assert (folder / f"{job.label}.xyz").is_file()
        assert (folder / "constraints.inp").is_file()
        assert (folder / f"{job.label}.out").is_file()

        inp = (folder / "constraints.inp").read_text()
        assert "$constrain" in inp
        assert "distance: 1, 2" in inp
        assert "force constant=0.5" in inp

        out = (folder / f"{job.label}.out").read_text()
        assert "--cinp constraints.inp" in out
        assert "--gfn2" in out
        assert "--alpb water" in out
